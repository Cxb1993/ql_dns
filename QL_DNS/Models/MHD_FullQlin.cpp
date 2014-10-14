//
//  MHD_FullQlin.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 9/16/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "MHD_FullQlin.h"

// Includes interaction of fluctuations with By and that's all
MHD_FullQlin::MHD_FullQlin(const Inputs& sp, MPIdata& mpi, fftwPlans& fft) :
equations_name("MHD_FullQlin"),
numMF_(4), numLin_(4),
q_(sp.q),
nu_(sp.nu), eta_(sp.eta),
f_noise_(sp.f_noise),QL_YN_(sp.QuasiLinearQ),
dont_drive_ky0_modes_(0),// If true, no driving ky=0 modes
Model(sp.NZ, sp.NXY , sp.L), // Dimensions - stored in base
mpi_(mpi), // MPI data
fft_(fft) // FFT data
{
    // Check that model is that specified in input file
    if (sp.equations_to_use != equations_name) {
        std::stringstream error_str;
        error_str << "Model name, " << equations_name << ", does not match that specified in input file, " << sp.equations_to_use << "!!" << std::endl;
        mpi.print1( error_str.str() );
        ABORT;
    }
    
    // Setup MPI
    mpi_.Split_NXY_Grid( Dimxy_full() ); // Work out MPI splitting
    // Assign K data
    K = new Kdata(this, &mpi_); // Stores all the K data
    // Fourier transform plans
    fft_.calculatePlans( NZfull(), NZ() );
    
    
    // Random generator
    mt_ = boost::random::mt19937( static_cast<unsigned int>(clock()+mpi_.my_n_v()) );  // Seed is from time
    ndist_ = boost::random::normal_distribution<double>(0,f_noise_/sqrt(2)); // Normal distribution, standard deviation f_noise_ (factor sqrt(2) is since it's complex)
    noise_buff_len_ = NZ()-1;
    noise_buff_ = dcmplxVec(num_Lin()*noise_buff_len_);// NZ()-1 since ky=kz=0 mode is not driven
    
    // Sizes of various arrays used for normalizing things
    totalN2_ = N(0)*N(1)*N(2); // Total number of grid points
    totalN2_ = totalN2_*totalN2_; // Squared
    mult_noise_fac_ = 1.0/(16*32*32); // Defined so that consistent with (16, 32, 32) results
    mult_noise_fac_ = mult_noise_fac_*mult_noise_fac_;
    reynolds_stress_ft_fac = 1.0/(N(1)*N(0))/(N(1)*N(0));
    
    
    ///////////////////////////////////////
    //  NB: MANY OF THE TEMP ASSIGNMENTS IN THIS VERSION
    //      ARE QUITE MEMORY INEFFICIENT. IF IT TURNS OUT
    //      TO BE PROBLEMATIC, COULD REASONABLY IMPROVE TO
    //      USE LESS THAN HALF OF WHAT IT DOES CURRENTLY,
    //      SACRIFICING SOME CLARITY
    ////////////////////////////////////////
    
    
    // Reynolds stresses
    Bx_drive_ = dcmplxVec::Zero(NZ());
    By_drive_ = dcmplxVec::Zero(NZ());
    Ux_drive_ = dcmplxVec::Zero(NZ());
    Uy_drive_ = dcmplxVec::Zero(NZ());
    reynolds_z_tmp_ = dcmplxVec::Zero(NZfull());
    // Send/receive buffers for MPI, for some reason MPI::IN_PLACE is not giving the correct result!
    num_to_mpi_send_ = (NZ()-1)/2+1; // Send only half of dealiased, fourier transformed reynolds stress (since it must be symmetric)
    reynolds_stress_MPI_send_ = dcmplxVec::Zero(num_MFs()*num_to_mpi_send_);
    reynolds_stress_MPI_receive_ = dcmplxVec::Zero(num_MFs()*num_to_mpi_send_);
    
    
    // CFL calculation
    kmax = 0;
    double dealiasfac[3] = {2.0, 2.0,3.0};
    for (int i=0; i<3; ++i) {
        if (kmax < 2*PI/L(i)*(N(i)/dealiasfac[i]))
            kmax = 2*PI/L(i)*(N(i)/dealiasfac[i]);
    }
    
    
    ////////////////////////////////////////////////////
    //               TEMPORARY ARRAYS                 //
    //  These are used during evaluation of the various parts of the model
    // There should never be any need to keep their value over more than 1 time-step
    //////////////////////////////////////////////////
    // These arrays are used for all sets of equations
    lapFtmp_ = doubVec( NZ() ); //For (time-dependent) k^2
    ilapFtmp_ = doubVec( NZ() ); // Inverse
    lap2tmp_ = doubVec( NZ() ); // For ky^2+kz^2 - could be pre-assigned
    ilap2tmp_ = doubVec( NZ() ); // Inverse
    
    // Temps for evaulation of main equations
    // u and b and derivatives
    u_ = new dcmplxVec[9];
    b_ = new dcmplxVec[9];
    for (int i=0; i<9; ++i) {
        u_[i] = dcmplxVec( NZfull() );
        b_[i] = dcmplxVec( NZfull() );
    }
    // Size NZfull()
    tmp1_z_ = dcmplxVec( NZfull() );
    // Size NZ() (kz,ky,kz)
    tmp1_k_ = dcmplxVec( NZ() );
    tmp2_k_ = dcmplxVec( NZ() );
    tmp3_k_ = dcmplxVec( NZ() );
    //  Size NZ() uy, uz etc.
    fluct_y_ = dcmplxVec( NZ() ); // k space versions of flutctuating y and z
    fluct_z_ = dcmplxVec( NZ() );



    
    //  Real space versions of mean fields - non-dealiased dimension
    MBx_ = dcmplxVec::Zero( NZfull() );
    dzMBx_ = dcmplxVec::Zero( NZfull() );
    MBy_ = dcmplxVec::Zero( NZfull() );
    dzMBy_ = dcmplxVec::Zero( NZfull() );
    MUx_ = dcmplxVec::Zero( NZfull() );
    dzMUx_ = dcmplxVec::Zero( NZfull() );
    MUy_ = dcmplxVec::Zero( NZfull() );
    dzMUy_ = dcmplxVec::Zero( NZfull() );

    
}


MHD_FullQlin::~MHD_FullQlin(){
    delete K;
    delete[] u_;
    delete[] b_;
}

void MHD_FullQlin::rhs(const double t, const double dt_lin,
                          const solution * SolIn, solution * SolOut,doubVec **linOpFluct) {
    // Calculate mean fields in real space
    // Calculate MFs in real space     By_ = MFin[1].matrix()*fft1Dfac_; // fft back doesn't include normalization
    // (NB could be optimized slightly by including memory copy in multiplication)
    
    // MEAN FIELDS IN REAL SPACE (Use tmp1_k_ as tmp variable)
    // Mean fields stored as (Bx, By, Ux, Uy)
    // Bx
    tmp1_k_ = (*SolIn->pMF(0))*fft_.fac1D();
    fft_.inverse( &tmp1_k_, &MBx_);
    tmp1_k_ = fft_.fac1D()*K->kz*(*SolIn->pMF(0));
    fft_.inverse( &tmp1_k_, &dzMBx_);
    // By
    tmp1_k_ = (*SolIn->pMF(1))*fft_.fac1D();
    fft_.inverse( &tmp1_k_, &MBy_);
    tmp1_k_ = fft_.fac1D()*K->kz*(*SolIn->pMF(1));
    fft_.inverse( &tmp1_k_, &dzMBy_);
    // Ux
    tmp1_k_ = (*SolIn->pMF(2))*fft_.fac1D();
    fft_.inverse( &tmp1_k_, &MUx_);
    tmp1_k_ = fft_.fac1D()*K->kz*(*SolIn->pMF(2));
    fft_.inverse( &tmp1_k_, &dzMUx_);
    // Uy
    tmp1_k_ = (*SolIn->pMF(3))*fft_.fac1D();
    fft_.inverse( &tmp1_k_, &MUy_);
    tmp1_k_ = fft_.fac1D()*K->kz*(*SolIn->pMF(3));
    fft_.inverse( &tmp1_k_, &dzMUy_);

//    //////////////
//    //  TO DELETE
//    MBx_.setZero();
//    dzMBx_.setZero();
//    MUx_.setZero();
//    dzMUx_.setZero();
//    MUy_.setZero();
//    dzMUy_.setZero();
//    //////////////

    
    // Reynolds stresses -- added to at each step in the loop
    reynolds_stress_MPI_send_.setZero();

    
    /////////////////////////////////////
    //   ALL OF THIS IS PARALLELIZED
    for (int i=0; i<Dimxy(); ++i) {
        // Full Loop containing all the main work equation
        ///////////////////////////////////////
        ///// MAIN LINEAR EQUATIONS
        
        assign_laplacians_(i, t, 1); // Assign kx, lapF etc.
        
        // u_ and b_ arrays in Real Space:
        // 0->ux, 1->uy, 2->uz, 3->dx(ux), 4->dx(uy), 5->dx(uz), 6->dy(ux), 7->dy(uy), 8->dy(uz)
        
        // Define uy, uz etc.,
        fluct_y_ = ( (-kyctmp_*kxctmp_)*(*SolIn->pLin(i, 0)) +  K->kz*(*SolIn->pLin(i, 1)) )*ilap2tmp_;
        fluct_z_ = ( -kxctmp_*K->kz*(*SolIn->pLin(i, 0)) -  kyctmp_*(*SolIn->pLin(i, 1)) )*ilap2tmp_;
        
        // Define derivatives and take ffts
        // u
        tmp1_k_ = fft_.fac1D()*(*SolIn->pLin(i, 0)); // ux (only for reynolds stress)
        fft_.inverse( &tmp1_k_, u_+0);
        tmp1_k_ = fft_.fac1D()*fluct_y_; // uy - (only for reynolds stress)
        fft_.inverse( &tmp1_k_, u_+1);
        tmp1_k_ = fft_.fac1D()*fluct_z_; // uz
        fft_.inverse( &tmp1_k_, u_+2);
        tmp1_k_ = (fft_.fac1D()*kxctmp_)*(*SolIn->pLin(i, 0)); // dx(ux)
        fft_.inverse( &tmp1_k_, u_+3); // dx(ux)
        tmp1_k_ = (fft_.fac1D()*kxctmp_)*fluct_y_;// dx(uy)
        fft_.inverse( &tmp1_k_, u_+4);
        tmp1_k_ = (fft_.fac1D()*kxctmp_)*fluct_z_;// dx(uz)
        fft_.inverse( &tmp1_k_, u_+5);
        tmp1_k_ = (fft_.fac1D()*kyctmp_)*(*SolIn->pLin(i, 0)); // dy(ux)
        fft_.inverse( &tmp1_k_, u_+6);
        tmp1_k_ = (fft_.fac1D()*kyctmp_)*fluct_y_;// dy(uy)
        fft_.inverse( &tmp1_k_, u_+7);
        tmp1_k_ = (fft_.fac1D()*kyctmp_)*fluct_z_;// dy(uz)
        fft_.inverse( &tmp1_k_, u_+8);
        
        // Define by, bz etc.,
        fluct_y_ = ( (-kyctmp_*kxctmp_)*(*SolIn->pLin(i, 2)) +  K->kz*(*SolIn->pLin(i, 3)) )*ilap2tmp_;
        fluct_z_ = ( -kxctmp_*K->kz*(*SolIn->pLin(i, 2)) -  kyctmp_*(*SolIn->pLin(i, 3)) )*ilap2tmp_;
        // b
        tmp1_k_ = fft_.fac1D()*(*SolIn->pLin(i, 2)); // bx (only for reynolds stress)
        fft_.inverse( &tmp1_k_, b_+0);
        tmp1_k_ = fft_.fac1D()*fluct_y_; // by - (only for reynolds stress)
        fft_.inverse( &tmp1_k_, b_+1);
        tmp1_k_ = fft_.fac1D()*fluct_z_; // bz
        fft_.inverse( &tmp1_k_, b_+2);
        tmp1_k_ = (fft_.fac1D()*kxctmp_)*(*SolIn->pLin(i, 2)); // dx(bx)
        fft_.inverse( &tmp1_k_, b_+3);
        tmp1_k_ = (fft_.fac1D()*kxctmp_)*fluct_y_;// dx(by)
        fft_.inverse( &tmp1_k_, b_+4);
        tmp1_k_ = (fft_.fac1D()*kxctmp_)*fluct_z_;// dx(bz)
        fft_.inverse( &tmp1_k_, b_+5);
        tmp1_k_ = (fft_.fac1D()*kyctmp_)*(*SolIn->pLin(i, 2)); // dy(bx)
        fft_.inverse( &tmp1_k_, b_+6);
        tmp1_k_ = (fft_.fac1D()*kyctmp_)*fluct_y_;// dy(by)
        fft_.inverse( &tmp1_k_, b_+7);
        tmp1_k_ = (fft_.fac1D()*kyctmp_)*fluct_z_;// dy(bz)
        fft_.inverse( &tmp1_k_, b_+8);
        
        ///////////////////////////////////
        //  Advance u and zeta
        // Form u.grad(u)-b.grad(b) - put into (tmp1_k_,tmp2_k_,tmp3_k_)
        // Ux dx(ux)+ Uy dy(ux) + uz dzUx
        tmp1_z_ = MUx_*u_[3] + MUy_*u_[6] + u_[2]*dzMUx_ -
                    (MBx_*b_[3] + MBy_*b_[6] + b_[2]*dzMBx_);
        fft_.forward(&tmp1_z_, &tmp1_k_);
        // Ux dx(uy)+ Uy dy(uy) + uz dzUy
        tmp1_z_ = MUx_*u_[4] + MUy_*u_[7] + u_[2]*dzMUy_ -
                    (MBx_*b_[4] + MBy_*b_[7] + b_[2]*dzMBy_);
        fft_.forward(&tmp1_z_, &tmp2_k_);
        // Ux dx(uz) + Uy dy(uz)
        tmp1_z_ = MUx_*u_[5] + MUy_*u_[8] - (MBx_*b_[5] + MBy_*b_[8]);
        fft_.forward(&tmp1_z_, &tmp3_k_);
        // u and zeta dot
        // u
        *SolOut->pLin(i, 0) = (  (-2.0*kxctmp_*kyctmp_*q_)*(*SolIn->pLin(i, 0)) +  2*K->kz*(*SolIn->pLin(i, 1))  +  (kxctmp_*kxctmp_)*tmp1_k_ + (kxctmp_*kyctmp_)*tmp2_k_ + kxctmp_*K->kz*tmp3_k_ )*ilapFtmp_    -     tmp1_k_;
        // zeta
        *SolOut->pLin(i, 1) = (q_-2.0)*K->kz*(*SolIn->pLin(i, 0))  -   (K->kz*tmp2_k_ - kyctmp_*tmp3_k_);
        ////////////////////////////////
        ///////////////////////////////////
        //  Advance b and eta
        // Form u.grad(b)-b.grad(u)   - put into (tmp1_k_,tmp2_k_,tmp3_k_)
        // Ux dx(bx)+ Uy dy(bx) + uz dzBx
        tmp1_z_ = MUx_*b_[3] + MUy_*b_[6] + u_[2]*dzMBx_ -
                    (MBx_*u_[3] + MBy_*u_[6] + b_[2]*dzMUx_);
        fft_.forward(&tmp1_z_, &tmp1_k_);
        // Ux dx(by)+ Uy dy(by) + uz dzBy
        tmp1_z_ = MUx_*b_[4] + MUy_*b_[7] + u_[2]*dzMBy_ -
                    (MBx_*u_[4] + MBy_*u_[7] + b_[2]*dzMUy_);
        fft_.forward(&tmp1_z_, &tmp2_k_);
        // Ux dx(bz) + Uy dy(bz)
        tmp1_z_ = MUx_*b_[5] + MUy_*b_[8] - (MBx_*u_[5] + MBy_*u_[8]);
        fft_.forward(&tmp1_z_, &tmp3_k_);
        // b and eta dot
        // b
        *SolOut->pLin(i, 2) = -tmp1_k_;
        // eta
        *SolOut->pLin(i, 3) = -q_*K->kz*(*SolIn->pLin(i, 2))   -   (K->kz*tmp2_k_ - kyctmp_*tmp3_k_);
        
        
        ///////////////////////////////////////
        ///// REYNOLDS STRESS
        // Chagned a little from other versions (MRIDSS)
        // Here, calculate all FFTs before doing MPI all reduce (and also multiply by necessary factors and kz to take derivative). This allows only half of the dealiased vector to be sent (at the expense of more ffts, reducing the communication load by 2/3
        
        // mult_fac for sum - since no -ve ky values are used, count everything but ky=0 twice
        double mult_fac = 2.0;
        if (kytmp_== 0.0 )
            mult_fac = 1.0; // Only count ky=0 mode once
        
        double ftfac = mult_fac*reynolds_stress_ft_fac; // Factor for summing all modes (1/(nx*ny)^2
        
        // Bx
        // bz*ux-uz*bx
        reynolds_z_tmp_ = ((b_[2]*u_[0].conjugate()).real()  -   (u_[2]*b_[0].conjugate()).real()).cast<dcmplx>();
        // Back to Fourier space
        fft_.forward(&reynolds_z_tmp_, &tmp1_k_);
        // Add on to MPI send vector
        reynolds_stress_MPI_send_.segment(0, num_to_mpi_send_) += ftfac*(K->kz*tmp1_k_).segment(0, num_to_mpi_send_);
        
        // By
        // bz*uy-uz*by
        reynolds_z_tmp_ = ((b_[2]*u_[1].conjugate()).real()  -   (u_[2]*b_[1].conjugate()).real()).cast<dcmplx>();
        // Back to fourier space
        fft_.forward(&reynolds_z_tmp_, &tmp1_k_);
        // Add on to MPI send vector
        reynolds_stress_MPI_send_.segment(num_to_mpi_send_, num_to_mpi_send_) += ftfac*(K->kz*tmp1_k_).segment(0, num_to_mpi_send_);
        
        // Ux
        // -(ux dx(ux) + uy dy(ux) + uz dz(ux)) + (bx dx(bx) + by dy(bx) + bz dz(bx))
        tmp1_k_ = fft_.fac1D()*K->kz*(*SolIn->pLin(i, 0)); // Need dz(ux)
        tmp2_k_ = fft_.fac1D()*K->kz*(*SolIn->pLin(i, 2)); // Need dz(bx)
        fft_.inverse(&tmp1_k_, &tmp1_z_); // Transforms
        fft_.inverse(&tmp2_k_, &reynolds_z_tmp_);
        reynolds_z_tmp_ = ( -((u_[0]*u_[3].conjugate()).real() + (u_[1]*u_[6].conjugate()).real() + (u_[2]*tmp1_z_.conjugate()).real() )   +
                ((b_[0]*b_[3].conjugate()).real() + (b_[1]*b_[6].conjugate()).real() + (b_[2]*reynolds_z_tmp_.conjugate()).real() )  ).cast<dcmplx>();
        // Back to Fourier space
        fft_.forward(&reynolds_z_tmp_, &tmp1_k_);
        // Add on to MPI send vector
        reynolds_stress_MPI_send_.segment(2*num_to_mpi_send_, num_to_mpi_send_) += ftfac*tmp1_k_.segment(0, num_to_mpi_send_);
        
        
        // Uy
        // -(ux dx(uy) + uy dy(uy) + uz dz(uy)) + (bx dx(by) + by dy(by) + bz dz(by))
        fluct_z_ = ( (-kyctmp_*kxctmp_)*(*SolIn->pLin(i, 0)) +  K->kz*(*SolIn->pLin(i, 1)) )*ilap2tmp_; // Redefine fluct_z to uy - fluct_y is already by
        tmp1_k_ = fft_.fac1D()*K->kz*fluct_z_; // Need dz(uy)
        tmp2_k_ = fft_.fac1D()*K->kz*fluct_y_; // Need dz(by)
        fft_.inverse(&tmp1_k_, &tmp1_z_); // Transforms
        fft_.inverse(&tmp2_k_, &reynolds_z_tmp_);
        reynolds_z_tmp_ = ( -((u_[0]*u_[4].conjugate()).real() + (u_[1]*u_[7].conjugate()).real() + (u_[2]*tmp1_z_.conjugate()).real() )   +
                        ((b_[0]*b_[4].conjugate()).real() + (b_[1]*b_[7].conjugate()).real() + (b_[2]*reynolds_z_tmp_.conjugate()).real() )  ).cast<dcmplx>();
        // Back to real space
        fft_.forward(&reynolds_z_tmp_, &tmp1_k_);
        // Add on to MPI send vector
        reynolds_stress_MPI_send_.segment(3*num_to_mpi_send_, num_to_mpi_send_) += ftfac*tmp1_k_.segment(0, num_to_mpi_send_);
        
        
        ///// END - REYNOLDS STRESS
        ///////////////////////////////////////
        
        ////////////////////////////////////////
        //////   LINEAR PART
        // Need to re-evaluate laplacian, since different time.
        kxtmp_=kxtmp_ + q_*dt_lin*kytmp_;
        linOpFluct[i][0] = nu_*((-kxtmp_*kxtmp_-kytmp_*kytmp_)+ K->kz2);
        linOpFluct[i][1] = linOpFluct[i][0];
        linOpFluct[i][2]  = (eta_/nu_)*linOpFluct[i][0];// Saves some calculation
        linOpFluct[i][3] = linOpFluct[i][2];
        
        ////////////////////////////////////////
        
        
    }
    //////////////////////////////////////
    // Reynolds stress
    // Sum accross all processes

    mpi_.SumAllReduce_dcmplx(reynolds_stress_MPI_send_.data(), reynolds_stress_MPI_receive_.data(), num_MFs()*num_to_mpi_send_);
    // Put into variables for MF update - flip to ensure realty
    // Bx
    Bx_drive_.segment(0, num_to_mpi_send_) = reynolds_stress_MPI_receive_.segment(0, num_to_mpi_send_);
    Bx_drive_.segment(num_to_mpi_send_,num_to_mpi_send_-1) = reynolds_stress_MPI_receive_.segment(1, num_to_mpi_send_-1).reverse().conjugate();
    // By
    By_drive_.segment(0, num_to_mpi_send_) = reynolds_stress_MPI_receive_.segment(num_to_mpi_send_, num_to_mpi_send_);
    By_drive_.segment(num_to_mpi_send_,num_to_mpi_send_-1) = reynolds_stress_MPI_receive_.segment(num_to_mpi_send_+1, num_to_mpi_send_-1).reverse().conjugate();
    // Ux
    Ux_drive_.segment(0, num_to_mpi_send_) = reynolds_stress_MPI_receive_.segment(2*num_to_mpi_send_, num_to_mpi_send_);
    Ux_drive_.segment(num_to_mpi_send_,num_to_mpi_send_-1) = reynolds_stress_MPI_receive_.segment(2*num_to_mpi_send_+1, num_to_mpi_send_-1).reverse().conjugate();
    // Uy
    Uy_drive_.segment(0, num_to_mpi_send_) = reynolds_stress_MPI_receive_.segment(3*num_to_mpi_send_, num_to_mpi_send_);
    Uy_drive_.segment(num_to_mpi_send_,num_to_mpi_send_-1) = reynolds_stress_MPI_receive_.segment(3*num_to_mpi_send_+1, num_to_mpi_send_-1).reverse().conjugate();
    
    // Reynolds stress
    //////////////////////////////////////
    
    
//    Bx_drive_.setZero();
//    By_drive_.setZero();
//    Ux_drive_.setZero();
//    Uy_drive_.setZero();
    //////////////////////////////////////
    ////   MEAN FIELDS    ////////////////
    if (QL_YN_) {
        *SolOut->pMF(1) = -q_*(*SolIn->pMF(0)) + By_drive_;
        *SolOut->pMF(0) = Bx_drive_;
        *SolOut->pMF(2) = 2*(*SolIn->pMF(3)) + Ux_drive_;
        *SolOut->pMF(3) = (q_-2.0)*(*SolIn->pMF(2)) + Uy_drive_;
    } else { // Still calculate Reynolds stress in linear calculation
        SolOut->pMF(1)->setZero();
        SolOut->pMF(0)->setZero();
        SolOut->pMF(2)->setZero();
        SolOut->pMF(3)->setZero();
    }

    
    
}

void MHD_FullQlin::linearOPs_Init(double t0, doubVec **linOpFluct, doubVec *linOpMF){
    
    // Fluctuations
    for (int i=0; i<Dimxy(); ++i) {
        assign_laplacians_(i, t0, 0); // Assign kx, lapF etc.
        linOpFluct[i][0] = nu_*lapFtmp_;
        linOpFluct[i][1] = linOpFluct[i][0];
        linOpFluct[i][2]  = (eta_/nu_)*linOpFluct[i][0];// Saves some calculation
        linOpFluct[i][3] = linOpFluct[i][2];
    }
    // Mean fields
    if (QL_YN_) {
        linOpMF[0]=eta_*K->kz2;
        linOpMF[1]=eta_*K->kz2;
        linOpMF[2]=nu_*K->kz2;
        linOpMF[3]=nu_*K->kz2;
    } else {
        linOpMF[0].setZero();
        linOpMF[1].setZero();
        linOpMF[2].setZero();
        linOpMF[3].setZero();
    }
    

}



//////////////////////////////////////////////////////////
//   Remapping
// Remapping for the shearing box
// Using the Lithwick method of continuously remapping wavenumbers when they become too large
void MHD_FullQlin::ShearingBox_Remap(double qt, solution *sol){
    double kxLfac=2*PI/L(0);// Define these for clarity
    int nx = N(0)-1;
    
    double kxt;
    for (int i=0; i<Dimxy(); ++i) {
        kxt = K->kx[i].imag() + qt*K->ky[i].imag();
        // If kx has grown too large
        if (kxt > nx/2.0*kxLfac) {
            // Put kx back in correct range
            K->kx[i] = K->kx[i] - dcmplx(0,nx*kxLfac);
            // zero out solution
            for (int Vn=0; Vn<num_Lin(); ++Vn) {
                sol->pLin(i, Vn)->setZero();
            }
        }
    }
}
//////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////
//    ENERGY/MOMENTUM ETC.
void MHD_FullQlin::Calc_Energy_AM_Diss(TimeVariables* tv, double t, const solution *sol) {
    // Energy, angular momentum and dissipation of the solution MFin and Cin
    // TimeVariables class stores info and data about energy, AM etc.
    // t is time
    
    // OUTPUT: energy[1] and [2] contain U and B mean field energies (energy[1]=0 for this)
    // energy[3] and [4] contain u and b fluctuating energies.
    
    
    tv->start_timing();
    //////////////////////////////////////////
    //    Can add more timevar outputs here
    //    Add to parallel loop below and update num_to_mpi
    
    // Energy
    double energy_u=0, energy_b=0;
    // Angular momentum
    double AM_u = 0, AM_b = 0;
    // Dissipation
    double diss_u=0,diss_b=0;
    
    // MPI buffers - this method passes around some unecessary data but will be very minimal
    const int num_to_mpi = 6;
    
    
    int mult_fac;// Factor to account for only doing half fft sum.
    /////////////////////////////////////
    //   ALL OF THIS IS PARALLELIZED
    for (int i=0; i<Dimxy();  ++i){
        // Form Laplacians using time-dependent kx
        assign_laplacians_(i, t, 1);
        
        mult_fac = 2;
        if (kytmp_== 0.0 )
            mult_fac = 1; // Only count ky=0 mode once
        
        lap2tmp_ = lapFtmp_*ilap2tmp_; // lap2tmp_ just a convenient storage
        
        if (tv->energy_save_Q()){
            //////////////////////////////////////
            //     ENERGY
            
            // (NB: ilap2tmp_ is negative but want abs, just use -ilap2)
            energy_u += mult_fac*(  lap2tmp_*( *(sol->pLin(i,0)) ).abs2() - ilap2tmp_*( *(sol->pLin(i,1)) ).abs2()  ).sum();
            energy_b += mult_fac*(  lap2tmp_*( *(sol->pLin(i,2)) ).abs2() - ilap2tmp_*( *(sol->pLin(i,3)) ).abs2()  ).sum();
            //////////////////////////////////////
        }
        if (tv->AngMom_save_Q()){
            //////////////////////////////////////
            //   Angular momentum
            
            fluct_y_ = ( (-kyctmp_*kxctmp_)*(*sol->pLin(i, 0)) +  K->kz*(*sol->pLin(i, 1)) )*ilap2tmp_;  //   Actually uy
            fluct_z_ = ( (-kyctmp_*kxctmp_)*(*sol->pLin(i, 2)) +  K->kz*(*sol->pLin(i, 3)) )*ilap2tmp_;  //   Actually by
            
            AM_u += mult_fac*( (*sol->pLin(i, 0))*fluct_y_.conjugate() ).real().sum();
            AM_b += mult_fac*( (*sol->pLin(i, 2))*fluct_z_.conjugate() ).real().sum();
            //
            //////////////////////////////////////
        }
        if (tv->dissip_save_Q()){
            //////////////////////////////////////
            //     DISSIPATION
            // (NB: ilap2tmp_ and lapFtmp_ are negative but want abs, just use negative)
            diss_u += (mult_fac*nu_)*(  -lapFtmp_*( lap2tmp_*( *(sol->pLin(i,0)) ).abs2() - ilap2tmp_*( *(sol->pLin(i,1)) ).abs2() ) ).sum();
            diss_b += (mult_fac*eta_)*(  -lapFtmp_*( lap2tmp_*( *(sol->pLin(i,2)) ).abs2() - ilap2tmp_*( *(sol->pLin(i,3)) ).abs2() ) ).sum();
            //
            //////////////////////////////////////
        }
        
        
    }
    double mpi_send_buff[num_to_mpi] = {energy_u,energy_b,AM_u,AM_b,diss_u,diss_b};
    double mpi_receive_buff[num_to_mpi];
    
    // Put the everything on processor 0
    mpi_.SumReduce_doub(mpi_send_buff,mpi_receive_buff,num_to_mpi);
    //    mpi_.SumReduce_IP_doub(&energy_u_f,1); // Is this working?
    
    // Currently, TimeVariables is set to save on root process, may want to generalize this at some point (SumReduce_doub is also)
    if (mpi_.my_n_v() == 0) {
        ////////////////////////////////////////////
        ///// All this is only on processor 0  /////
        double divfac=1.0/totalN2_;
        double divfavMF = 1.0/( NZfull()*NZfull() );
        energy_u = mpi_receive_buff[0]*divfac;
        energy_b = mpi_receive_buff[1]*divfac;
        AM_u = mpi_receive_buff[2]*divfac;
        AM_b = mpi_receive_buff[3]*divfac;
        diss_u = mpi_receive_buff[4]*divfac;
        diss_b = mpi_receive_buff[5]*divfac;
        
        
        ///////////////////////////////////////
        ///       MEAN FIELDS            //////
        // Only need to calculate on one processor
        double energy_MU=0, energy_MB=0;
        double AM_MU =0, AM_MB=0;
        double diss_MU =0, diss_MB =0;
        
        energy_MB = ( (*(sol->pMF(0))).abs2().sum() + (*(sol->pMF(1))).abs2().sum() )*divfavMF;
        energy_MU = ( (*(sol->pMF(2))).abs2().sum() + (*(sol->pMF(3))).abs2().sum() )*divfavMF;
        
        AM_MB = divfavMF*( (*sol->pMF(0))*(sol->pMF(1)->conjugate()) ).real().sum();
        AM_MU = divfavMF*( (*sol->pMF(2))*(sol->pMF(3)->conjugate()) ).real().sum();
        
        diss_MB = (divfavMF*eta_)*( -(K->kz2*sol->pMF(0)->abs2()).sum() - (K->kz2*sol->pMF(1)->abs2()).sum() );
        diss_MU = (divfavMF*eta_)*( -(K->kz2*sol->pMF(2)->abs2()).sum() - (K->kz2*sol->pMF(3)->abs2()).sum() );

        
        ///////////////////////////////////////
        //////         OUTPUT            //////
        
        // Energy
        double* en_point = tv->current_energy();
        en_point[0] = energy_MU/2;
        en_point[1] = energy_MB/2;
        en_point[2] = energy_u/2;
        en_point[3] = energy_b/2;
        
        // Angular momentum
        double* AM_point = tv->current_AM();
        AM_point[0] = AM_MU;
        AM_point[1] = AM_MB;
        AM_point[2] = AM_u;
        AM_point[3] = AM_b;
        
        // Energy
        double* diss_point = tv->current_diss();
        diss_point[0] = diss_MU;
        diss_point[1] = diss_MB;
        diss_point[2] = diss_u;
        diss_point[3] = diss_b;
        
        ///// All this is only on processor 0  /////
        ////////////////////////////////////////////
    }
    
    tv->Save_Data(t);
    tv->finish_timing();

    
}


//////////////////////////////////////////////////////////
//   AUXILIARY FUNCTIONS
inline void MHD_FullQlin::assign_laplacians_(int i, double t, bool need_inverse){
    ind_ky_ = K->ky_index[i];
    // Form Laplacians using time-dependent kx
    kyctmp_ = K->ky[i];
    kytmp_ = kyctmp_.imag(); // It is useful to have both complex and double versions
    kxctmp_ = K->kx[i] + q_*t*kyctmp_;
    kxtmp_ = kxctmp_.imag();
    
    lap2tmp_ = K->lap2[ind_ky_];
    lapFtmp_ = -kxtmp_*kxtmp_ + lap2tmp_;
    if (need_inverse){
        ilap2tmp_ = K->ilap2[ind_ky_];
        ilapFtmp_ = 1/lapFtmp_;
        if (ind_ky_ == 0 && K->kx_index[i]==0)
            ilapFtmp_(0) = 1.0; // Otherwise get nans
    }
}
//////////////////////////////////////////////////////////



//////////////////////////
// CFL number
double MHD_FullQlin::Calculate_CFL(const solution *sol)  {
// Returns CFL/dt to calculate dt - in this case CFL/dt = kmax By + q
    
    // Mean fields have been previously calculated, so may as well use them
    double Bxmax = sqrt(MBx_.abs2().maxCoeff());
    double Bymax = sqrt(MBy_.abs2().maxCoeff());
    double Uxmax = sqrt(MUx_.abs2().maxCoeff());
    double Uymax = sqrt(MUy_.abs2().maxCoeff());
    // CFL
    return kmax*(Bymax + Bxmax + Uxmax + Uymax) + q_;
    
}



//////////////////////////////////////////////////////////
//   DRIVING NOISE
void MHD_FullQlin::DrivingNoise(double t, double dt, solution *sol) {
    // NZ() is now dealiased NZ!!!
    // So - drive all of NZ and Nyquist frequency is not included
    
    // Adds noise onto linear part of solution
    for (int i=0; i<Dimxy(); ++i) {
        
        if (K->ky_index[i] != 0){ // ky=0 dealt with seperately. kx=ky=0 missed automatically
            
            ///////////////////////////
            // Noise multipliers
            assign_laplacians_(i,t,0); // Assign laplacian's - no inverses
            
            // ky=0 so not many if statements
            double noise_multfac = dt*totalN2_*mult_noise_fac_; // f_noise is included in normal distribution
            lapFtmp_ = noise_multfac*lap2tmp_/lapFtmp_; // lapFtmp_ is no longer lapF!
            
            lap2tmp_ = (-noise_multfac*lap2tmp_).sqrt();
            lapFtmp_ = lapFtmp_.sqrt();
            
            double *multU_pnt = lapFtmp_.data();
            double *multZeta_pnt = lap2tmp_.data();
            /////////////////////////////
            
            // ky != 0 modes are completely random
            // U
            dcmplx* ePoint = (*(sol->pLin(i,0) )).data();
            for (int jj=0; jj<NZ(); ++jj) {
                ePoint[jj] += multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
            }
            // zeta
            ePoint = (*(sol->pLin(i,1) )).data();
            for (int jj=0; jj<NZ(); ++jj) {
                ePoint[jj] += multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
            }
            // b
            ePoint = (*(sol->pLin(i,2) )).data();
            for (int jj=0; jj<NZ(); ++jj) {
                ePoint[jj] += multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
            }
            // eta
            ePoint = (*(sol->pLin(i,3) )).data();
            for (int jj=0; jj<NZ(); ++jj) {
                ePoint[jj] += multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
            }
        }
        
        // Don't drive kx=ky=0 modes (mean fields) or ky=kz=0 modes (non-invertible).These are both taken care of in laplacian
        // Still need to make sure that ky=0 kx!=0 kz!=0 modes are symmetric. In particular ky=0, kx=a, kz=b must be the complex conjugate to kx=-a, kz=-b
        // To do this, have to communicate between processors!
        
        
    }
    
    if (!dont_drive_ky0_modes_){
        // Deal with ky=0 part - a real PIA!!
        // Iterators for vectors
        K->i_tosend = K->match_kx_tosend.begin();
        K->i_loc = K->match_kx_local.begin();
        K->i_from = K->match_kx_from.begin();
        K->i_fp = K->match_kx_from_p.begin();
        K->i_sp = K->match_kx_tosend_p.begin();
        
#ifdef USE_MPI_FLAG
        //////////////////////////////////////////////////
        // Sending - create noise here and broadcast
        while (K->i_sp < K->match_kx_tosend_p.end()) {
            
            int i = *(K->i_tosend); // Index in k array
            if (K->ky_index[i]!=0) // Be paranoid
                mpi_.print1("Warning: Something wrong in noise data swapping!!");
            ///////////////////////////
            // Noise multipliers
            
            assign_laplacians_(i,t,0); // Assign laplacian's - no inverses
            
            double noise_multfac = dt*totalN2_*mult_noise_fac_; // f_noise is included in normal distribution
            lap2tmp_(0) = 0.0; // Don't drive ky=kz=0 (not really necessary)
            
            lapFtmp_ = noise_multfac*lap2tmp_/lapFtmp_; // lapFtmp_ is no longer lapF!
            
            lap2tmp_ = (-noise_multfac*lap2tmp_).sqrt();
            lapFtmp_ = lapFtmp_.sqrt();
            
            double *multU_pnt = lapFtmp_.data();
            double *multZeta_pnt = lap2tmp_.data();
            /////////////////////////////
            
            // Make some noise - add to current k value
            dcmplx* noise_buff_dat_ = noise_buff_.data();
            dcmplx* ePoint = (*(sol->pLin(i,0) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                noise_buff_dat_[jj-1]=multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                ePoint[jj] += noise_buff_dat_[jj-1];
            }
            // zeta
            noise_buff_dat_ += noise_buff_len_;
            ePoint = (*(sol->pLin(i,1) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                noise_buff_dat_[jj-1]=multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                ePoint[jj] += noise_buff_dat_[jj-1];
            }
            // b
            noise_buff_dat_ += noise_buff_len_;
            ePoint = (*(sol->pLin(i,2) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                noise_buff_dat_[jj-1]=multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                ePoint[jj] += noise_buff_dat_[jj-1];
            }
            // eta
            noise_buff_dat_ += noise_buff_len_;
            ePoint = (*(sol->pLin(i,3) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                noise_buff_dat_[jj-1]=multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                ePoint[jj] += noise_buff_dat_[jj-1];
            }
            
            // Send data
            mpi_.Send_dcmplx(noise_buff_.data(), num_Lin()*noise_buff_len_, *(K->i_sp), *(K->i_tosend));
            
            //        cout << "(kx,ky)=(" << K->kx[i] <<"," << K->ky[i] << ")" << endl;
            //        cout << noise_buff_.segment(0, nz).transpose() << endl<<endl;
            // Update iterator
            ++(K->i_tosend);
            ++(K->i_sp);
        }
        ////////////////////////////////////////////
        // Recieving
        while (K->i_fp < K->match_kx_from_p.end()) {
            
            int i = *(K->i_loc);
            // Receive data from matching call
            mpi_.Recv_dcmplx(noise_buff_.data(), num_Lin()*noise_buff_len_, *(K->i_fp), *(K->i_from));
            
            // Flip noise around
            for (int nV=0; nV<num_Lin(); ++nV) {
                noise_buff_.segment(nV*noise_buff_len_, noise_buff_len_).reverseInPlace();
            }
            
            //        cout << "(kx,ky)=(" << K->kx[i] <<"," << K->ky[i] << ")" << endl;
            //        cout << noise_buff_.segment(0, nz).transpose().conjugate() <<endl<<endl;
            
            // Add to solution
            for (int nV=0; nV<num_Lin(); ++nV) {
                (*(sol->pLin(i,nV) )).segment(1,NZ()-1) += noise_buff_.segment(nV*noise_buff_len_, noise_buff_len_).conjugate();
            }
            
            
            // Update iterators
            ++(K->i_fp);
            ++(K->i_loc);
            ++(K->i_from);
        }
        
        // mpi_.Barrier();  // This is probably not necessary but don't see why it would cause harm
        ////////////////////////////////////////////////////////////////////
#else       // Need a separate "no-mpi" case, since MPI_Send functions not defined
        ////////////////////////////////////////////////////////////////////
        
        while (K->i_tosend < K->match_kx_tosend.end()) {
            ///////////////////////////
            // Noise multipliers
            int i = *(K->i_tosend); // index
            assign_laplacians_(i,t,0); // Assign laplacian's - no inverses
            
            double noise_multfac = dt*totalN2_*mult_noise_fac_; // f_noise is included in normal distribution
            lap2tmp_(0) = 0.0; // Don't drive ky=kz=0
            
            lapFtmp_ = noise_multfac*lap2tmp_/lapFtmp_; // lapFtmp_ is no longer lapF!
            
            lap2tmp_ = (-noise_multfac*lap2tmp_).sqrt();
            lapFtmp_ = lapFtmp_.sqrt();
            
            double *multU_pnt = lapFtmp_.data();
            double *multZeta_pnt = lap2tmp_.data();
            /////////////////////////////
            
            // Make some noise - add to current k value
            dcmplx* noise_buff_dat_ = noise_buff_.data();
            dcmplx* ePoint = (*(sol->pLin(i,0) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                noise_buff_dat_[jj-1]=multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                ePoint[jj] += noise_buff_dat_[jj-1];
            }
            // zeta
            noise_buff_dat_ += noise_buff_len_;
            ePoint = (*(sol->pLin(i,1) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                noise_buff_dat_[jj-1]=multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                ePoint[jj] += noise_buff_dat_[jj-1];
            }
            // b
            noise_buff_dat_ += noise_buff_len_;
            ePoint = (*(sol->pLin(i,2) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                noise_buff_dat_[jj-1]=multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                ePoint[jj] += noise_buff_dat_[jj-1];
            }
            // eta
            noise_buff_dat_ += noise_buff_len_;
            ePoint = (*(sol->pLin(i,3) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                noise_buff_dat_[jj-1]=multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                ePoint[jj] += noise_buff_dat_[jj-1];
            }
            
            
            // Find matching K value and use the same noise
            i = *(K->i_loc);
            
            // Flip noise around
            for (int nV=0; nV<num_Lin(); ++nV) {
                noise_buff_.segment(nV*noise_buff_len_, noise_buff_len_).reverseInPlace();
            }
            
            // Add to solution
            for (int nV=0; nV<num_Lin(); ++nV) {
                (*(sol->pLin(i,nV) )).segment(1,NZ()-1) += noise_buff_.segment(nV*noise_buff_len_, noise_buff_len_).conjugate();
            }
            
            
            // Update iterators
            ++(K->i_tosend);
            ++(K->i_loc);
        }
        
#endif
        ////////////////////////////////////////////////////////////////////
        
    }
    
}
//////////////////////////////////////////////////////////

