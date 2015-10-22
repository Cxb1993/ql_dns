//
//  HD_fullU.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 9/16/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "HD_fullU.h"

// Includes interaction of fluctuations with By and that's all
HD_fullU::HD_fullU(const Inputs& sp, MPIdata& mpi, fftwPlans& fft) :
equations_name("HD_fullU"),
numMF_(2), numLin_(2),
q_(sp.q),Omega_(sp.omega), // Rotation
nu_(sp.nu),
f_noise_(sp.f_noise),QL_YN_(sp.QuasiLinearQ),
dont_drive_ky0_modes_(0),// If true, no driving ky=0 modes
save_full_reynolds_(1),// Saves all the nonlinear stresses if true
Model(sp.NZ, sp.NXY , sp.L), // Dimensions - stored in base
mpi_(mpi), // MPI data
fft_(fft) // FFT data
{
    // Check that model is that specified in input file - unecessary now
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
    // Noise cutoff - compare to laplacian so squared
    noise_range_[0] = sp.noise_range_low*sp.noise_range_low;
    noise_range_[1] = sp.noise_range_high*sp.noise_range_high;
    drive_condition_ = Eigen::Matrix<bool,Eigen::Dynamic,1>(NZ());
    print_noise_range_();

    if (save_full_reynolds_)
        num_reynolds_saves_ = 2*NZfull();// ONLY U STRESSES
    else num_reynolds_saves_=5;
    
        
    
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
    
    // Random settings, print to avoid mistakes
    std::stringstream prnt;
    prnt << "Rotation set to " << Omega_ <<"\n";
    mpi_.print1(prnt.str());
    if (dont_drive_ky0_modes_)
        mpi_.print1("ky=0 modes are excluded from this calculation!\n");
    
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
    for (int i=0; i<9; ++i) {
        u_[i] = dcmplxVec( NZfull() );
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
    MUx_ = dcmplxVec::Zero( NZfull() );
    dzMUx_ = dcmplxVec::Zero( NZfull() );
    MUy_ = dcmplxVec::Zero( NZfull() );
    dzMUy_ = dcmplxVec::Zero( NZfull() );

    
}


HD_fullU::~HD_fullU(){
    delete K;
    delete[] u_;
}

void HD_fullU::rhs(const double t, const double dt_lin,
                          const solution * SolIn, solution * SolOut,doubVec **linOpFluct) {
    // Calculate mean fields in real space
    // Calculate MFs in real space     By_ = MFin[1].matrix()*fft1Dfac_; // fft back doesn't include normalization
    // (NB could be optimized slightly by including memory copy in multiplication)
    
    // MEAN FIELDS IN REAL SPACE (Use tmp1_k_ as tmp variable)
    // Mean fields stored as (Bx, By, Ux, Uy)
    // Ux
    tmp1_k_ = (*SolIn->pMF(0))*fft_.fac1D();
    fft_.inverse( &tmp1_k_, &MUx_);
    tmp1_k_ = fft_.fac1D()*K->kz*(*SolIn->pMF(0));
    fft_.inverse( &tmp1_k_, &dzMUx_);
    // Uy
    tmp1_k_ = (*SolIn->pMF(1))*fft_.fac1D();
    fft_.inverse( &tmp1_k_, &MUy_);
    tmp1_k_ = fft_.fac1D()*K->kz*(*SolIn->pMF(1));
    fft_.inverse( &tmp1_k_, &dzMUy_);



    
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
        
        
        ///////////////////////////////////
        //  Advance u and zeta
        // Form u.grad(u)-b.grad(b) - put into (tmp1_k_,tmp2_k_,tmp3_k_)
        // Ux dx(ux)+ Uy dy(ux) + uz dzUx
        tmp1_z_ = MUx_*u_[3] + MUy_*u_[6] + u_[2]*dzMUx_;
        fft_.forward(&tmp1_z_, &tmp1_k_);
        // Ux dx(uy)+ Uy dy(uy) + uz dzUy
        tmp1_z_ = MUx_*u_[4] + MUy_*u_[7] + u_[2]*dzMUy_ ;
        fft_.forward(&tmp1_z_, &tmp2_k_);
        // Ux dx(uz) + Uy dy(uz)
        tmp1_z_ = MUx_*u_[5] + MUy_*u_[8];
        fft_.forward(&tmp1_z_, &tmp3_k_);
        // u and zeta dot
        // u
        *SolOut->pLin(i, 0) = (  (-2.0*kxctmp_*kyctmp_*q_)*(*SolIn->pLin(i, 0)) +  2*Omega_*K->kz*(*SolIn->pLin(i, 1))  +  (kxctmp_*kxctmp_)*tmp1_k_ + (kxctmp_*kyctmp_)*tmp2_k_ + kxctmp_*K->kz*tmp3_k_ )*ilapFtmp_     -     tmp1_k_;
        // zeta
        *SolOut->pLin(i, 1) = (q_-2.0*Omega_)*K->kz*(*SolIn->pLin(i, 0))  -   (K->kz*tmp2_k_ - kyctmp_*tmp3_k_);
        ////////////////////////////////
        
        
        ///////////////////////////////////////
        ///// REYNOLDS STRESS
        // Changed a little from other versions (MRIDSS)
        // Here, calculate all FFTs before doing MPI all reduce (and also multiply by necessary factors and kz to take derivative). This allows only half of the dealiased vector to be sent (at the expense of more ffts, reducing the communication load by 2/3
        
        // mult_fac for sum - since no -ve ky values are used, count everything but ky=0 twice
        double mult_fac = 2.0;
        if (kytmp_== 0.0 )
            mult_fac = 1.0; // Only count ky=0 mode once
        
        double ftfac = mult_fac*reynolds_stress_ft_fac; // Factor for summing all modes (1/(nx*ny)^2
        
        
        // Ux
        // -(ux dx(ux) + uy dy(ux) + uz dz(ux))
        tmp1_k_ = fft_.fac1D()*K->kz*(*SolIn->pLin(i, 0)); // Need dz(ux)
        fft_.inverse(&tmp1_k_, &tmp1_z_); // Transforms
        reynolds_z_tmp_ = ( -((u_[0]*u_[3].conjugate()).real() + (u_[1]*u_[6].conjugate()).real() + (u_[2]*tmp1_z_.conjugate()).real() )   ).cast<dcmplx>();
        // Back to Fourier space
        fft_.forward(&reynolds_z_tmp_, &tmp1_k_);
        // Add on to MPI send vector
        reynolds_stress_MPI_send_.segment(0, num_to_mpi_send_) += ftfac*tmp1_k_.segment(0, num_to_mpi_send_);
        
        
        // Uy
        // -(ux dx(uy) + uy dy(uy) + uz dz(uy))
        fluct_z_ = ( (-kyctmp_*kxctmp_)*(*SolIn->pLin(i, 0)) +  K->kz*(*SolIn->pLin(i, 1)) )*ilap2tmp_; // Redefine fluct_z to uy - fluct_y is already by
        tmp1_k_ = fft_.fac1D()*K->kz*fluct_z_; // Need dz(uy)
        fft_.inverse(&tmp1_k_, &tmp1_z_); // Transforms
        reynolds_z_tmp_ = ( -((u_[0]*u_[4].conjugate()).real() + (u_[1]*u_[7].conjugate()).real() + (u_[2]*tmp1_z_.conjugate()).real() )  ).cast<dcmplx>();
        // Back to real space
        fft_.forward(&reynolds_z_tmp_, &tmp1_k_);
        // Add on to MPI send vector
        reynolds_stress_MPI_send_.segment(num_to_mpi_send_, num_to_mpi_send_) += ftfac*tmp1_k_.segment(0, num_to_mpi_send_);
        
        
        ///// END - REYNOLDS STRESS
        ///////////////////////////////////////
        
        ////////////////////////////////////////
        //////   LINEAR PART
        // Need to re-evaluate laplacian, since different time.
        kxtmp_=kxtmp_ + q_*dt_lin*kytmp_;
        linOpFluct[i][0] = nu_*((-kxtmp_*kxtmp_-kytmp_*kytmp_)+ K->kz2);
        linOpFluct[i][1] = linOpFluct[i][0];
        
        ////////////////////////////////////////
        
        
    }
    //////////////////////////////////////
    // Reynolds stress
    // Sum accross all processes

    mpi_.SumAllReduce_dcmplx(reynolds_stress_MPI_send_.data(), reynolds_stress_MPI_receive_.data(), num_MFs()*num_to_mpi_send_);
    // Put into variables for MF update - flip to ensure realty
    // Ux
    Ux_drive_.segment(0, num_to_mpi_send_) = reynolds_stress_MPI_receive_.segment(0, num_to_mpi_send_);
    Ux_drive_.segment(num_to_mpi_send_,num_to_mpi_send_-1) = reynolds_stress_MPI_receive_.segment(1, num_to_mpi_send_-1).reverse().conjugate();
    // Uy
    Uy_drive_.segment(0, num_to_mpi_send_) = reynolds_stress_MPI_receive_.segment(num_to_mpi_send_, num_to_mpi_send_);
    Uy_drive_.segment(num_to_mpi_send_,num_to_mpi_send_-1) = reynolds_stress_MPI_receive_.segment(num_to_mpi_send_+1, num_to_mpi_send_-1).reverse().conjugate();
    
    // Reynolds stress
    //////////////////////////////////////
    

    //////////////////////////////////////
    ////   MEAN FIELDS    ////////////////
    if (QL_YN_) {
        // U update
        *SolOut->pMF(0) = 2*Omega_*(*SolIn->pMF(1)) + Ux_drive_;
        *SolOut->pMF(1) = (q_-2.0*Omega_)*(*SolIn->pMF(0)) + Uy_drive_;
    } else { // Still calculate Reynolds stress in linear calculation
        SolOut->pMF(1)->setZero();
        SolOut->pMF(0)->setZero();
    }

    
    
}

void HD_fullU::linearOPs_Init(double t0, doubVec **linOpFluct, doubVec *linOpMF){
    
    // Fluctuations
    for (int i=0; i<Dimxy(); ++i) {
        assign_laplacians_(i, t0, 0); // Assign kx, lapF etc.
        linOpFluct[i][0] = nu_*lapFtmp_;
        linOpFluct[i][1] = linOpFluct[i][0];
    }
    // Mean fields
    if (QL_YN_) {
        linOpMF[0]=nu_*K->kz2;
        linOpMF[1]=nu_*K->kz2;
    } else {
        linOpMF[0].setZero();
        linOpMF[1].setZero();
    }
    

}



//////////////////////////////////////////////////////////
//   Remapping
// Remapping for the shearing box
// Using the Lithwick method of continuously remapping wavenumbers when they become too large
void HD_fullU::ShearingBox_Remap(double qt, solution *sol){
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
void HD_fullU::Calc_Energy_AM_Diss(TimeVariables* tv, double t, const solution *sol) {
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
    double energy_u=0;
    // Angular momentum
    double AM_u = 0;
    // Dissipation
    double diss_u=0;
    
    // MPI buffers - this method passes around some unecessary data but will be very minimal
    const int num_to_mpi = 3;
    
    
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
            //////////////////////////////////////
            
//            std::cout << "kx = " << kxtmp_ << ", ky = " << kytmp_ << ": " <<mult_fac*(  lap2tmp_*( *(sol->pLin(i,0)) ).abs2() - ilap2tmp_*( *(sol->pLin(i,1)) ).abs2()  ).sum() <<std::endl;
        }
        if (tv->AngMom_save_Q()){
            //////////////////////////////////////
            //   Angular momentum
            
            fluct_y_ = ( (-kyctmp_*kxctmp_)*(*sol->pLin(i, 0)) +  K->kz*(*sol->pLin(i, 1)) )*ilap2tmp_;  //   Actually uy
            
            AM_u += mult_fac*( (*sol->pLin(i, 0))*fluct_y_.conjugate() ).real().sum();
            //
            //////////////////////////////////////
        }
        if (tv->dissip_save_Q()){
            //////////////////////////////////////
            //     DISSIPATION
            // (NB: ilap2tmp_ and lapFtmp_ are negative but want abs, just use negative)
            diss_u += (mult_fac*nu_)*(  -lapFtmp_*( lap2tmp_*( *(sol->pLin(i,0)) ).abs2() - ilap2tmp_*( *(sol->pLin(i,1)) ).abs2() ) ).sum();
            
            //
            //////////////////////////////////////
        }
        
        
    }
    double mpi_send_buff[num_to_mpi] = {energy_u,AM_u,diss_u};
    double mpi_receive_buff[num_to_mpi];
    
    // Put the everything on processor 0
    mpi_.SumReduce_doub(mpi_send_buff,mpi_receive_buff,num_to_mpi);
    //    mpi_.SumReduce_IP_doub(&energy_u_f,1); // Is this working?
    

    // Currently, TimeVariables is set to save on root process, may want to generalize this at some point (SumReduce_doub is also)
    if (mpi_.my_n_v() == 0) {
        ////////////////////////////////////////////
        ///// All this is only on processor 0  /////
        double divfac=L(2)/totalN2_;
        
        energy_u = mpi_receive_buff[0]*divfac;
        AM_u = mpi_receive_buff[1]*divfac;
        diss_u = mpi_receive_buff[2]*divfac;
        
        ///////////////////////////////////////
        ///       MEAN FIELDS            //////
        // Only need to calculate on one processor
        double energy_MU=0;
        double AM_MU =0;
        double diss_MU =0;
        double divfavMF = L(2)/( NZfull()*NZfull() );
        
        energy_MU = ( (*(sol->pMF(0))).abs2().sum() + (*(sol->pMF(1))).abs2().sum() )*divfavMF;
        
        AM_MU = divfavMF*( (*sol->pMF(0))*(sol->pMF(1)->conjugate()) ).real().sum();
        
        diss_MU = (divfavMF*nu_)*( -(K->kz2*sol->pMF(0)->abs2()).sum() - (K->kz2*sol->pMF(1)->abs2()).sum() );
        
        
        ///////////////////////////////////////
        //////         OUTPUT            //////
        
        // Energy
        double* en_point = tv->current_energy();
        en_point[0] = energy_MU/2;
        en_point[1] = energy_u/2;
        
        // Angular momentum
        double* AM_point = tv->current_AM();
        AM_point[0] = AM_MU;
        AM_point[1] = AM_u;
        
        // Energy
        double* diss_point = tv->current_diss();
        diss_point[0] = diss_MU;
        diss_point[1] = diss_u;
        
        ///// All this is only on processor 0  /////
        ////////////////////////////////////////////
        
        // Set etaK_times_kmax in TimeVariables - This is then used in output for monitoring spatial resolution as simulation progresses
        // Use kz_max for this scale (sometimes have high x resolution)
        // Calculate on zero and broadcast
        tv->etaK_times_kmax[0] = pow(nu_*nu_*nu_/(diss_MU+diss_u),0.25)*(2*PI/L(2)*(N(2)/3));
        
        
        if (tv->reynolds_save_Q()) {
            //////////////////////////////////////
            //   Reynolds stress
            // Saves quantities to do with the reynolds stress and dynamo
            // 1) Shear contribution: Re( -q Bx By)/|By|
            // 2) y emf: Re( bzuy_m_uzby_c_*By )/|By|
            // 3) y dissipation: eta*k^2*By
            // 4) x emf: Re( bzux_m_uzbx_c_*Bx )/|Bx|
            // 5) x dissipation: eta*k^2*Bx
            // Have assumed k0 to be lowest kz! i.e., driving largest dynamo possible in the box
            // There may be slight errors here from saving using the updated values of MFin, presumably this is a small effect, especially at high resolution (low dt) and in steady state.
            double* rey_point = tv->current_reynolds();
            
            if (!save_full_reynolds_){
                // FOR MEASURING DYNAMO IN FIXED B - k1 is wavelength of B
                // NB: Had a mistake before, was 1/sqrt(sol->pMF(1)->abs2().sum()) rather than 1/(sol->pMF(1)->abs2().sum()). This makes the answer depend on the chosen value of the imposed field.
                
                
                mpi_.print1("REYNOLDS SAVE NOT PROPERLY IMPLEMENTED!!!!!");
                // Alpha effect (mean zero) - alpha_yy*k1^2
                rey_point[0] = (Ux_drive_*K->kz*sol->pMF(1)->conjugate()).real().sum()/(sol->pMF(1)->abs2().sum());
                // Shear current - k1^2*eta_yx
                rey_point[1] = (Bx_drive_*sol->pMF(1)->conjugate()).real().sum()/(sol->pMF(1)->abs2().sum());
                // Turbulent By diffusion - -k1^2*eta_yy
                rey_point[2] = (By_drive_*sol->pMF(1)->conjugate()).real().sum()/(sol->pMF(1)->abs2().sum());
                // Bx nonzero
                // Off diagonal - k1^2*eta_xy
                rey_point[3] = (By_drive_*sol->pMF(0)->conjugate()).real().sum()/(sol->pMF(0)->abs2().sum());
                // Turbulent Bx diffusion - -k1^2*eta_xx
                rey_point[4] = (Bx_drive_*sol->pMF(0)->conjugate()).real().sum()/sol->pMF(0)->abs2().sum();
                
 
            } else {
                // SAVE FULL NONLINEAR STRESS DATA
                // Only saving B data for now
                fft_.inverse( &Ux_drive_, &reynolds_z_tmp_);
                for (int jj=0; jj<NZfull(); ++jj) {
                    rey_point[jj] = real(reynolds_z_tmp_(jj))*fft_.fac1D();
                }
                fft_.inverse( &Uy_drive_, &reynolds_z_tmp_);
                for (int jj=0; jj<NZfull(); ++jj) {
                    rey_point[jj+NZfull()] = real(reynolds_z_tmp_(jj))*fft_.fac1D();
                }
                
                
            }

            //////////////////////////////////////
            
            
        }
        
    }
    
    tv->Save_Data(t);
    tv->finish_timing();
    
    mpi_.BroadcastFromNode0_doub(tv->etaK_times_kmax, 2);
//    if ( tv->etaK_times_kmax[0]<0.5 ||  tv->etaK_times_kmax[1]<0.5)
//        mpi_.print1("Warning: eta_K is not well resolved by spatial grid!!!\n");
    
    
}


//////////////////////////////////////////////////////////
//   AUXILIARY FUNCTIONS
inline void HD_fullU::assign_laplacians_(int i, double t, bool need_inverse){
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



void HD_fullU::print_noise_range_(){
    // Print out total number of driven modes
    int tot_count = 0, driv_count = 0;
    int tot_count_all = 0, driv_count_all = 0;// MPI reduced versions
    for (int i=0; i<Dimxy(); ++i) {
        assign_laplacians_(i, 0.0, 0);
        if (dont_drive_ky0_modes_ && kytmp_ == 0.0) {
            drive_condition_.setZero();
        } else {
            drive_condition_ = lapFtmp_> -noise_range_[1] && lapFtmp_< -noise_range_[0];
        }
        tot_count += NZ();
        driv_count += drive_condition_.cast<int>().sum();
    }
    mpi_.SumReduce_int(&tot_count, &tot_count_all, 1);
    mpi_.SumReduce_int(&driv_count, &driv_count_all, 1);
    std::stringstream noise_output;
    noise_output << "Driving " << driv_count_all << " out of a total of " << tot_count_all << " modes (k = " << sqrt(noise_range_[0]) << " to " << sqrt(noise_range_[1]) << ")\n";
    mpi_.print1(noise_output.str());
    
}
//////////////////////////////////////////////////////////



//////////////////////////
// CFL number
double HD_fullU::Calculate_CFL(const solution *sol)  {
// Returns CFL/dt to calculate dt - in this case CFL/dt = kmax U + q
// Including the mean field equation also (treating u.Gu as a simple source to dtU=-qU+S) adds on
// max(u.Gu/U) -- The quotient makes for a problem, so use mean(k*u.Gu)/mean(k*U) to weight towards high K
    
    // Mean fields have been previously calculated, so may as well use them
    double Uxmax = sqrt(MUx_.abs2().maxCoeff());
    double Uymax = sqrt(MUy_.abs2().maxCoeff());
    // This only includes the linear part, so...
    // Add in part using Ux_drive and Uy_drive directly
    double MFcont = 0.;
    if (QL_YN_) {
        // Include an energy cuttoff for incuding this part -- only important for large mean fields - pretty ad-hoc, but actually seems to work relatively well....
        if (( MUx_.abs2().sum() + MUy_.abs2().sum() )*L(2) >1e-6) {
            MFcont = (Ux_drive_*K->kz2).abs().mean()/((*sol->pMF(0))*K->kz2).abs().mean()  +  (Uy_drive_*K->kz2).abs().mean()/((*sol->pMF(0))*K->kz2).abs().mean();
        }
    }
    
    // CFL
//    std::stringstream out;
//    out << "kmax*(Uxmax + Uymax), " << kmax*(Uxmax + Uymax) << ": MFcont, " << MFcont <<  "\n";
//    mpi_.print1(out.str());
    
    
    return kmax*(Uxmax + Uymax) + q_ + MFcont;
}



//////////////////////////////////////////////////////////
//   DRIVING NOISE
void HD_fullU::DrivingNoise(double t, double dt, solution *sol) {
    // NZ() is now dealiased NZ!!!
    // So - drive all of NZ and Nyquist frequency is not included
    
    // Scale the magnetic driving by the value of the mean-field - tangling of the mean-field (see Yousef et al)
    // 0 turns this off, non-zero is the proportionality
    double drive_mag_mult = 1.0; // Scaling magnetic noise
    
    // Adds noise onto linear part of solution
    for (int i=0; i<Dimxy(); ++i) {
        
        if (K->ky_index[i] != 0){ // ky=0 dealt with seperately. kx=ky=0 missed automatically
            
            ///////////////////////////
            // Noise multipliers
            assign_laplacians_(i,t,0); // Assign laplacian's - no inverses
            
            // ky=0 so not many if statements
            double noise_multfac = dt*totalN2_*mult_noise_fac_; // f_noise is included in normal distribution
            // If k is outside of noise_range_, don't drive.
            drive_condition_ = lapFtmp_> -noise_range_[1] && lapFtmp_< -noise_range_[0]; // lapFtmp is negative, so conditions backwards
            
            lapFtmp_ = noise_multfac*lap2tmp_/lapFtmp_; // lapFtmp_ is no longer lapF!
            
            lap2tmp_ = (-noise_multfac*lap2tmp_).sqrt();
            lapFtmp_ = lapFtmp_.sqrt();
            
            double *multU_pnt = lapFtmp_.data();
            double *multZeta_pnt = lap2tmp_.data();
            bool *drive_cond_pnt = drive_condition_.data();

            /////////////////////////////
            
            // ky != 0 modes are completely random
            // U
            dcmplx* ePoint = (*(sol->pLin(i,0) )).data();
            for (int jj=0; jj<NZ(); ++jj) {
                if (drive_cond_pnt[jj])
                    ePoint[jj] += multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
            }
            // zeta
            ePoint = (*(sol->pLin(i,1) )).data();
            for (int jj=0; jj<NZ(); ++jj) {
                if (drive_cond_pnt[jj])
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
            // If k is outside of noise_range_, don't drive.
            drive_condition_ = lapFtmp_> -noise_range_[1] && lapFtmp_< -noise_range_[0]; // lapFtmp is negative, so conditions backwards
            
            lapFtmp_ = noise_multfac*lap2tmp_/lapFtmp_; // lapFtmp_ is no longer lapF!
            
            lap2tmp_ = (-noise_multfac*lap2tmp_).sqrt();
            lapFtmp_ = lapFtmp_.sqrt();
            
            double *multU_pnt = lapFtmp_.data();
            double *multZeta_pnt = lap2tmp_.data();
            bool *drive_cond_pnt = drive_condition_.data();
            /////////////////////////////
            
            // Make some noise - add to current k value
            noise_buff_.setZero();
            dcmplx* noise_buff_dat_ = noise_buff_.data();
            dcmplx* ePoint = (*(sol->pLin(i,0) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                if (drive_cond_pnt[jj]){
                    noise_buff_dat_[jj-1]=multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                    ePoint[jj] += noise_buff_dat_[jj-1];
                }
            }
            // zeta
            noise_buff_dat_ += noise_buff_len_;
            ePoint = (*(sol->pLin(i,1) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                if (drive_cond_pnt[jj]){
                    noise_buff_dat_[jj-1]=multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                    ePoint[jj] += noise_buff_dat_[jj-1];
                }
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
            // If k is outside of noise_range_, don't drive.
            drive_condition_ = lapFtmp_> -noise_range_[1] && lapFtmp_< -noise_range_[0]; // lapFtmp is negative, so conditions backwards
            
            lapFtmp_ = noise_multfac*lap2tmp_/lapFtmp_; // lapFtmp_ is no longer lapF!
            
            lap2tmp_ = (-noise_multfac*lap2tmp_).sqrt();
            lapFtmp_ = lapFtmp_.sqrt();
            
            double *multU_pnt = lapFtmp_.data();
            double *multZeta_pnt = lap2tmp_.data();
            bool *drive_cond_pnt = drive_condition_.data();
            /////////////////////////////
            
            // Make some noise - add to current k value
            dcmplx* noise_buff_dat_ = noise_buff_.data();
            dcmplx* ePoint = (*(sol->pLin(i,0) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                if (drive_cond_pnt[jj] ){
                    noise_buff_dat_[jj-1]=multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                    ePoint[jj] += noise_buff_dat_[jj-1];
                }
            }
            // zeta
            noise_buff_dat_ += noise_buff_len_;
            ePoint = (*(sol->pLin(i,1) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                if (drive_cond_pnt[jj] ){
                    noise_buff_dat_[jj-1]=multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                    ePoint[jj] += noise_buff_dat_[jj-1];
                }
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

