//
//  CompMHDIso_FixedMF.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 3/2/16.
//  Copyright © 2016 J Squire. All rights reserved.
//

#include "CompMHDIso_FixedMF.h"

// Compressible MHD
// Mean fields are rho, Ux, Uy, Bx, By (with constant Bz also), and full compressible equations in linear variables.
// NB: No MF stresses, fully linear
CompMHDIso_FixedMF::CompMHDIso_FixedMF(const Inputs& sp, MPIdata& mpi, fftwPlans& fft) :
equations_name("CompMHDIso_FixedMF"),
numMF_(4), numLin_(7),
q_(sp.q),Omega_(sp.omega), // Rotation
B0z_(sp.B0z), // Mean vertical field
rho0_(sp.rho0), // Mean density
P0_(sp.P0), // Background pressure, determines sound speed
cs20_(P0_/rho0_), // Sound speed
nu_(sp.nu), eta_(sp.eta), kappa_(sp.kappa),
viscosity_order_(sp.viscosity_order),
f_noise_(sp.f_noise),QL_YN_(sp.QuasiLinearQ),
dont_drive_ky0_modes_(0),// If true, no driving ky=0 modes
drive_only_velocity_fluctuations_(sp.drive_only_velocityQ), // If true, no magnetic forcing
drive_only_magnetic_fluctuations_(sp.drive_only_magneticQ), // If true, no velocity forcing
save_full_reynolds_(0),// Saves all the nonlinear stresses if true
apply_b_divergence_cleaning_(1), // Whether to clean div(b)
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
    
    // Linear variables are
    // 0: rho, 1: vx, 2: vy, 3: vz, 4: bx, 5: by, 6: bz
    
    // Setup MPI
    mpi_.Split_NXY_Grid( Dimxy_full() ); // Work out MPI splitting
    // Assign K data
    K = new Kdata(this, &mpi_); // Stores all the K data
    // Fourier transform plans
    fft_.calculatePlans( NZfull(), NZ() );
    
    
    
    
    // Random generator //static_cast<unsigned int>(clock()+mpi_.my_n_v())
    mt_ = boost::random::mt19937( static_cast<unsigned int>(clock()+mpi_.my_n_v()) );  // Seed is from time
    ndist_ = boost::random::normal_distribution<double>(0,f_noise_/sqrt(2)); // Normal distribution, standard deviation f_noise_ (factor sqrt(2) is since it's complex)
    noise_buff_len_ = NZ()-1;// NZ()-1 since ky=kz=0 mode is not driven
    num_noise_vars_ = 4;
    noise_buff_ = dcmplxVec(num_noise_vars_*noise_buff_len_); // Generate divergence free noise so only 4 cmpts
    // Noise cutoff - compare to laplacian so squared
    noise_range_[0] = sp.noise_range_low*sp.noise_range_low;
    noise_range_[1] = sp.noise_range_high*sp.noise_range_high;
    drive_condition_ = Eigen::Matrix<bool,Eigen::Dynamic,1>(NZ());
    print_noise_range_();
    
    if (save_full_reynolds_)
        num_reynolds_saves_ = 4*NZfull();
    else num_reynolds_saves_=4;   // Saves stress on each MF
    
    
    
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
    
    // Random settings, print to avoid mistakes
    std::stringstream prnt;
    prnt << "Omega = " << Omega_ << ", shear = " << q_ << "\n";
    mpi_.print1(prnt.str());
    if (dont_drive_ky0_modes_)
        mpi_.print1("ky=0 modes are excluded from this calculation!\n");
    if (drive_only_velocity_fluctuations_)
        mpi_.print1("No Magnetic driving noise!\n");
    if (drive_only_magnetic_fluctuations_)
        mpi_.print1("No velocity driving noise!\n");
    if (!apply_b_divergence_cleaning_)
        mpi_.print1("No magnetic field divergence cleaning!\n");
    // Viscosity at kmax, useful for hyperviscosities
    prnt.str("");
    prnt << "nu*k^n ~ 1 at k = " << pow(1./nu_,1.0/viscosity_order_) << " (" << pow(1./nu_,1./viscosity_order_)/kmax << "*kmax): nu*kmax^n = " << nu_*pow(kmax,viscosity_order_) << std::endl << std::endl;
    mpi_.print1(prnt.str());
    
    
    
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
        u_[i].setZero();
        b_[i] = dcmplxVec( NZfull() );
        b_[i].setZero();
    }
    rho_ = dcmplxVec( NZfull() );
    divu_ = dcmplxVec( NZfull() );
    // Size NZfull()
    tmp1_z_ = dcmplxVec( NZfull() );
    // Size NZ() (kz,ky,kz)
    tmp1_k_ = dcmplxVec( NZ() );
    tmp2_k_ = dcmplxVec( NZ() );
    tmp3_k_ = dcmplxVec( NZ() );
    tmp4_k_ = dcmplxVec( NZ() );
    

    
    
    
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


CompMHDIso_FixedMF::~CompMHDIso_FixedMF(){
    delete K;
    delete[] u_;
    delete[] b_;
}

// Modified const of SolIn so that it can be divergence cleaned!

void CompMHDIso_FixedMF::rhs(const double t, const double dt_lin,
                        solution * SolIn, solution * SolOut,doubVec **linOpFluct) {
    // Calculate mean fields in real space
    // Calculate MFs in real space     By_ = MFin[1].matrix()*fft1Dfac_; // fft back doesn't include normalization
    // (NB could be optimized slightly by including memory copy in multiplication)
    
    // MEAN FIELDS IN REAL SPACE (Use tmp1_k_ as tmp variable)
    // Mean fields stored as (Bx, By, Ux, Uy)
    // Even though these are fixed in time, keep the sam format for passing things around here.
    // Assumed rho0 is constant
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
    
    
    
    
    // Reynolds stresses -- added to at each step in the loop
    reynolds_stress_MPI_send_.setZero();
    divb_monitor_ = 0.;
    divb_global_ = 0.;// to sum for mpi
    
    /////////////////////////////////////
    //   ALL OF THIS IS PARALLELIZED
    for (int i=0; i<Dimxy(); ++i) {
        // Full Loop containing all the main work equation
        ///////////////////////////////////////
        ///// MAIN LINEAR EQUATIONS
        
        assign_laplacians_(i, t, 1); // Assign kx, lapF etc.
        
        // u_ and b_ arrays in Real Space:
        // 0->ux, 1->uy, 2->uz, 3->dx(ux), 4->dx(uy), 5->dx(uz), 6->dy(ux), 7->dy(uy), 8->dy(uz), 9->div(u)
        
        // Linear variables are
        // 0: rho, 1: vx, 2: vy, 3: vz, 4: bx, 5: jx
        
        
        // Define uy, uz etc.
        
        // Define derivatives and take ffts -- whole bunch of these are a bit pointless, since don't need to take fft to take x or y deriv... should change later!
        // u
        tmp1_k_ = fft_.fac1D()*(*SolIn->pLin(i, 1)); // ux (only for reynolds stress)
        fft_.inverse( &tmp1_k_, u_+0);
        tmp1_k_ = fft_.fac1D()*(*SolIn->pLin(i, 2)); // uy - (only for reynolds stress)
        fft_.inverse( &tmp1_k_, u_+1);
        tmp1_k_ = fft_.fac1D()*(*SolIn->pLin(i, 3)); // uz
        fft_.inverse( &tmp1_k_, u_+2);
        u_[3] = kxctmp_*u_[0];// dx(ux)
        u_[4] = kxctmp_*u_[1];// dx(uy)
        u_[5] = kxctmp_*u_[2];// dx(uz)
        u_[6] = kyctmp_*u_[0];// dy(ux)
        u_[7] = kyctmp_*u_[1];// dy(uy)
        u_[8] = kyctmp_*u_[2];// dy(uz)
        // And div(u) in real space
        tmp1_k_ = fft_.fac1D()*(kxctmp_*(*SolIn->pLin(i, 1)) + kyctmp_*(*SolIn->pLin(i, 2)) + K->kz*(*SolIn->pLin(i, 3)));
        fft_.inverse( &tmp1_k_, &divu_);
//        cout << divu_.transpose().abs().sum() << endl;

        
        // Clean divergence from b -- put into (*SolIn->pLin(i, 4)), (*SolIn->pLin(i, 5)) etc.. Not the most efficient, but easy and obvious and meh
        // Without shear, this is conserved. With shear, only to time integrator precision, presumably because kx changes each time step. This then quickly causes divergence to grow because in the formulation I have here, dt(divb)!=0 if input divb!=0 (would need to also include U divB term, then it wouldn't grow so fast).
        if (apply_b_divergence_cleaning_) {
            tmp1_k_ = ilapFtmp_*(kxctmp_*(*SolIn->pLin(i, 4)) + kyctmp_*(*SolIn->pLin(i, 5)) + K->kz*(*SolIn->pLin(i, 6)));
            (*SolIn->pLin(i, 4)) -= kxctmp_*tmp1_k_;
            (*SolIn->pLin(i, 5)) -= kyctmp_*tmp1_k_;
            (*SolIn->pLin(i, 6)) -= K->kz*tmp1_k_;
        }
        // Monitor div(b) to make sure it doens't grow too large overall
        divb_monitor_ += (kxctmp_*(*SolIn->pLin(i, 4))+kyctmp_*(*SolIn->pLin(i, 5))+K->kz*(*SolIn->pLin(i, 6))).abs().sum()/  max(1e-6,sqrt(SolIn->pLin(i, 4)->abs2().sum()+SolIn->pLin(i, 5)->abs2().sum()+SolIn->pLin(i, 6)->abs2().sum()));

        
        
        // b
        tmp1_k_ = fft_.fac1D()*(*SolIn->pLin(i, 4)); // bx (only for reynolds stress)
        fft_.inverse( &tmp1_k_, b_+0);
        tmp1_k_ = fft_.fac1D()*(*SolIn->pLin(i, 5)); // by - (only for reynolds stress)
        fft_.inverse( &tmp1_k_, b_+1);
        tmp1_k_ = fft_.fac1D()*(*SolIn->pLin(i, 6)); // bz
        fft_.inverse( &tmp1_k_, b_+2);
        b_[3] = kxctmp_*b_[0];
        b_[4] = kxctmp_*b_[1];
        b_[5] = kxctmp_*b_[2];
        b_[6] = kyctmp_*b_[0];
        b_[7] = kyctmp_*b_[1];
        b_[8] = kyctmp_*b_[2];

        
        // Linear variables are
        // 0: rho, 1: vx, 2: vy, 3: vz, 4: bx, 5: by, 6: bz
        
        ///////////////////////////////////
        //  Advance rho
        tmp1_k_ = fft_.fac1D()*(*SolIn->pLin(i, 0));
        fft_.inverse( &tmp1_k_, &rho_); //  Real space rho
        // U.Grho
        tmp1_z_ = (MUx_*kxctmp_ + MUy_*kyctmp_)*rho_;
        fft_.forward(&tmp1_z_, &tmp1_k_);
        *SolOut->pLin(i, 0) = -tmp1_k_ - rho0_*(kxctmp_*(*SolIn->pLin(i, 1)) + kyctmp_*(*SolIn->pLin(i, 2)) + K->kz*(*SolIn->pLin(i, 3)));
        
//        cout << "kx,ky: " << kxtmp_ << ", " << kytmp_ << " t=" << t << endl;
//        cout << rho_.real().transpose() << endl;
        
        ///////////////////////////////////
        //  Advance u
        // Form u.grad(u)-b.grad(b) - put into (tmp1_k_,tmp2_k_,tmp3_k_)
        // Ux dx(ux)+ Uy dy(ux) + uz dzUx - Bterms
        tmp1_z_ = MUx_*u_[3] + MUy_*u_[6] + u_[2]*dzMUx_ -  (MBx_*b_[3] + MBy_*b_[6] + b_[2]*dzMBx_ )/rho0_;
        fft_.forward(&tmp1_z_, &tmp1_k_);
        // Ux dx(uy)+ Uy dy(uy) + uz dzUy - Bterms
        tmp1_z_ = MUx_*u_[4] + MUy_*u_[7] + u_[2]*dzMUy_ - (MBx_*b_[4] + MBy_*b_[7] + b_[2]*dzMBy_ )/rho0_;
        fft_.forward(&tmp1_z_, &tmp2_k_);
        // Ux dx(uz) + Uy dy(uz)- Bterms
        tmp1_z_ = MUx_*u_[5] + MUy_*u_[8] - (MBx_*b_[5] + MBy_*b_[8])/rho0_;
        fft_.forward(&tmp1_z_, &tmp3_k_);
        // Magnetic pressure b.B - take gradient in equations.
        tmp1_z_ = MBx_*b_[0] + MBy_*b_[1] + B0z_*b_[2];
        fft_.forward(&tmp1_z_, &tmp4_k_);
        // u dot
        *SolOut->pLin(i, 1) = - tmp1_k_ + 2.0*Omega_*(*SolIn->pLin(i, 2)) + B0z_/rho0_*K->kz*(*SolIn->pLin(i, 4))  - 1/rho0_*kxctmp_*(cs20_*(*SolIn->pLin(i, 0)) + tmp4_k_);
        *SolOut->pLin(i, 2) = -tmp2_k_ + (q_-2.0*Omega_)*(*SolIn->pLin(i, 1)) + B0z_/rho0_*K->kz*(*SolIn->pLin(i, 5)) - 1/rho0_*kyctmp_*(cs20_*(*SolIn->pLin(i, 0)) + tmp4_k_);
        *SolOut->pLin(i, 3) = -tmp3_k_  + B0z_/rho0_*K->kz*(*SolIn->pLin(i, 6)) - 1/rho0_*K->kz*(cs20_*(*SolIn->pLin(i, 0)) + tmp4_k_);
        
//        cout << (rho0_*K->kz*(cs20_*(*SolIn->pLin(i, 0)))).transpose() << endl;

//        // As a test, clean divergence
//        tmp4_k_ = ilapFtmp_*(kxctmp_*(*SolOut->pLin(i, 1)) + kyctmp_*(*SolOut->pLin(i, 2)) + K->kz*(*SolOut->pLin(i, 3)));
//        (*SolOut->pLin(i, 1)) -=  kxctmp_*tmp4_k_;
//        (*SolOut->pLin(i, 2)) -=  kyctmp_*tmp4_k_;
//        (*SolOut->pLin(i, 3)) -=  K->kz*tmp4_k_;
        
        
        ////////////////////////////////
        // u_ and b_ arrays in Real Space:
        // 0->ux, 1->uy, 2->uz, 3->dx(ux), 4->dx(uy), 5->dx(uz), 6->dy(ux), 7->dy(uy), 8->dy(uz)
        ///////////////////////////////////
        //  Advance b and eta
        // Form u.grad(b)-b.grad(u)   - put into (tmp1_k_,tmp2_k_,tmp3_k_)
        // Ux dx(bx)+ Uy dy(bx) + uz dzBx - B.GU + B*div(u)
        tmp1_z_ = MUx_*b_[3] + MUy_*b_[6] + u_[2]*dzMBx_ - (MBx_*u_[3] + MBy_*u_[6] + b_[2]*dzMUx_) + MBx_*divu_;
        fft_.forward(&tmp1_z_, &tmp1_k_);
        // Ux dx(by)+ Uy dy(by) + uz dzBy - B.GU + B*div(u)
        tmp1_z_ = MUx_*b_[4] + MUy_*b_[7] + u_[2]*dzMBy_ - (MBx_*u_[4] + MBy_*u_[7] + b_[2]*dzMUy_) + MBy_*divu_;
        fft_.forward(&tmp1_z_, &tmp2_k_);
        // Ux dx(bz) + Uy dy(bz) - B.GU + B*div(u)
        tmp1_z_ = MUx_*b_[5] + MUy_*b_[8] - (MBx_*u_[5] + MBy_*u_[8])  +  B0z_*divu_;
        fft_.forward(&tmp1_z_, &tmp3_k_);
        // bx, by, bz dot
        // bx
        *SolOut->pLin(i, 4) = B0z_*K->kz*(*SolIn->pLin(i, 1)) - tmp1_k_;
        // by
        *SolOut->pLin(i, 5) = -q_*(*SolIn->pLin(i, 4))  +   B0z_*K->kz*(*SolIn->pLin(i, 2))  -   tmp2_k_;
        // bz
        *SolOut->pLin(i, 6) =  B0z_*K->kz*(*SolIn->pLin(i, 3))  -  tmp3_k_;


//        cout << "kx, ky:" << kxtmp_ <<" " << kytmp_ << ":\n";
//        cout << divu_.transpose() << endl;
//        cout << (kxctmp_*(*SolIn->pLin(i, 4))+kyctmp_*(*SolIn->pLin(i, 5))+K->kz*(*SolIn->pLin(i, 6))).transpose() << endl;
//        cout << (kxctmp_*(*SolOut->pLin(i, 4))+kyctmp_*(*SolOut->pLin(i, 5))+K->kz*(*SolOut->pLin(i, 6))).transpose() << endl <<endl;

        
       
        
        ////////////////////////////////////////
        //////   LINEAR PART
        // Need to re-evaluate laplacian, since different time.
        kxtmp_=kxtmp_ + q_*dt_lin/2*kytmp_;
        // Use lapFtmp_ to store Laplacian^Viscosity_order/2
        lapFtmp_= -((-kxtmp_*kxtmp_-kytmp_*kytmp_)+ K->kz2).abs().pow(viscosity_order_/2);
        linOpFluct[i][0] = kappa_*lapFtmp_;
        linOpFluct[i][1] = nu_*lapFtmp_;
        linOpFluct[i][2] = linOpFluct[i][1];
        linOpFluct[i][3] = linOpFluct[i][1];
        linOpFluct[i][4]  = eta_*lapFtmp_;
        linOpFluct[i][5] = linOpFluct[i][4];
        linOpFluct[i][6] = linOpFluct[i][4];
        
        ////////////////////////////////////////
        
        
        
        
        
    }
    
    mpi_.SumAllReduce_doub(&divb_monitor_,&divb_global_,1);
    if (divb_global_>1e-9){
        stringstream prnt;
        prnt << "div(b)/b has grown to " << divb_global_ << " at t="<< t <<". Be worried!\n";
        mpi_.print1(prnt.str());
    }
    
    // No Reynolds stress
    //////////////////////////////////////
    
    
    //////////////////////////////////////
    ////   MEAN FIELDS    ////////////////
    if (QL_YN_) {
        mpi_.print1("CompMHDIso_FixedMF is not set up for mean field evolution!!\n");
    } else { // Still calculate Reynolds stress in linear calculation
        SolOut->pMF(1)->setZero();
        SolOut->pMF(0)->setZero();
        SolOut->pMF(2)->setZero();
        SolOut->pMF(3)->setZero();
    }
    
    
    
}

void CompMHDIso_FixedMF::linearOPs_Init(double t0, doubVec **linOpFluct, doubVec *linOpMF){
    
    // Fluctuations
    for (int i=0; i<Dimxy(); ++i) {
        assign_laplacians_(i, t0, 0); // Assign kx, lapF etc.
        // Use lapFtmp_ to store Laplacian^Viscosity_order/2
        lapFtmp_= -lapFtmp_.abs().pow(viscosity_order_/2);
        linOpFluct[i][0] = kappa_*lapFtmp_;
        linOpFluct[i][1] = nu_*lapFtmp_;
        linOpFluct[i][2] = linOpFluct[i][1];
        linOpFluct[i][3] = linOpFluct[i][1];
        linOpFluct[i][4]  = eta_*lapFtmp_;
        linOpFluct[i][5] = linOpFluct[i][4];
        linOpFluct[i][6] = linOpFluct[i][4];

    }
    // Mean fields
    if (QL_YN_) {
        mpi_.print1("CompMHDIso_FixedMF is not set up for mean field evolution!!\n\n");
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
void CompMHDIso_FixedMF::ShearingBox_Remap(double qt, solution *sol){
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
void CompMHDIso_FixedMF::Calc_Energy_AM_Diss(TimeVariables* tv, double t, const solution *sol) {
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
        assign_laplacians_(i, t, 0);
        
        mult_fac = 2;
        if (kytmp_== 0.0 )
            mult_fac = 1; // Only count ky=0 mode once
        
        
        if (tv->energy_save_Q()){
            //////////////////////////////////////
            //     ENERGY
            
            energy_u += mult_fac*rho0_*( (*(sol->pLin(i,1))).abs2() + (*(sol->pLin(i,2))).abs2() + (*(sol->pLin(i,3))).abs2() ).sum();
            energy_b += mult_fac*( (*(sol->pLin(i,4))).abs2() + (*(sol->pLin(i,5))).abs2() + (*(sol->pLin(i,6))).abs2() ).sum();
            //////////////////////////////////////
        }
        if (tv->AngMom_save_Q()){
            //////////////////////////////////////
            //   Angular momentum
            
            AM_u += mult_fac*( (*(sol->pLin(i,1)))*(*(sol->pLin(i,2))).conjugate() ).real().sum();
            AM_b += mult_fac*( (*sol->pLin(i, 4))*(*sol->pLin(i, 5)).conjugate() ).real().sum();
            //
            //////////////////////////////////////
        }
        if (tv->dissip_save_Q()){
            //////////////////////////////////////
            //     DISSIPATION
            // (NB: ilap2tmp_ and lapFtmp_ are negative but want abs, just use negative)
            diss_u += (mult_fac*nu_)*(  -lapFtmp_*( (*(sol->pLin(i,1))).abs2() + (*(sol->pLin(i,2))).abs2() + (*(sol->pLin(i,3))).abs2()  ) ).sum();
            diss_b += (mult_fac*eta_)*(  -lapFtmp_*( (*(sol->pLin(i,4))).abs2() + (*(sol->pLin(i,5))).abs2() + (*(sol->pLin(i,6))).abs2() ) ).sum();
            //
            //////////////////////////////////////
        }
        
        if (tv->reynolds_save_Q()){
            ///////////////////////////////////////
            ///// REYNOLDS STRESS
            // No MF evolution, so only calculate these when it's actually necessary
            
            double ftfac = mult_fac*reynolds_stress_ft_fac; // Factor for summing all modes (1/(nx*ny)^2
            
            // u
            tmp1_k_ = fft_.fac1D()*(*sol->pLin(i, 1)); // ux
            fft_.inverse( &tmp1_k_, u_+0);
            tmp1_k_ = fft_.fac1D()*(*sol->pLin(i, 2)); // uy
            fft_.inverse( &tmp1_k_, u_+1);
            tmp1_k_ = fft_.fac1D()*(*sol->pLin(i, 3)); // uz
            fft_.inverse( &tmp1_k_, u_+2);
            // And div(u) in real space
            tmp1_k_ = fft_.fac1D()*(kxctmp_*(*sol->pLin(i, 1)) + kyctmp_*(*sol->pLin(i, 2)) + K->kz*(*sol->pLin(i, 3)));
            fft_.inverse( &tmp1_k_, &divu_);
            
            // b
            tmp1_k_ = fft_.fac1D()*(*sol->pLin(i, 4)); // bx
            fft_.inverse( &tmp1_k_, b_+0);
            tmp1_k_ = fft_.fac1D()*(*sol->pLin(i, 5)); // by
            fft_.inverse( &tmp1_k_, b_+1);
            tmp1_k_ = fft_.fac1D()*(*sol->pLin(i, 6)); // bz
            fft_.inverse( &tmp1_k_, b_+2);
            
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
            
            // Ux -- probably a better way of working these out...
            // -(ux dx(ux) + uy dy(ux) + uz dz(ux)) + (bx dx(bx) + by dy(bx) + bz dz(bx))
            // Could use -dz(uz ux)+(ux divu) + dz(bz bx)instead
            reynolds_z_tmp_ = (-u_[2]*u_[0].conjugate().real() + b_[2]*b_[0].conjugate().real()).cast<dcmplx>();
            tmp1_z_ = (divu_*u_[0].conjugate()).real().cast<dcmplx>();
            // Back to Fourier space
            fft_.forward(&reynolds_z_tmp_, &tmp1_k_);
            fft_.forward(&tmp1_z_, &tmp2_k_);
            // Add on to MPI send vector
            reynolds_stress_MPI_send_.segment(2*num_to_mpi_send_, num_to_mpi_send_) += ftfac*(K->kz*tmp1_k_ + tmp2_k_).segment(0, num_to_mpi_send_);
            
            
            // Uy
            // -(ux dx(uy) + uy dy(uy) + uz dz(uy)) + (bx dx(by) + by dy(by) + bz dz(by))
            // Could use -dz(uz uy)+(uy divu) + dz(bz by)instead
            reynolds_z_tmp_ = (-u_[2]*u_[1].conjugate().real() + b_[2]*b_[1].conjugate().real()).cast<dcmplx>();
            tmp1_z_ = (divu_*u_[1].conjugate()).real().cast<dcmplx>();
            // Back to Fourier space
            fft_.forward(&reynolds_z_tmp_, &tmp1_k_);
            fft_.forward(&tmp1_z_, &tmp2_k_);
            // Add on to MPI send vector
            reynolds_stress_MPI_send_.segment(3*num_to_mpi_send_, num_to_mpi_send_) += ftfac*(K->kz*tmp1_k_ + tmp2_k_).segment(0, num_to_mpi_send_);
            
            
            ///// END - REYNOLDS STRESS
            ///////////////////////////////////////
            
            
//            //////////////////
//            // TO DELETE!!!!
//            tmp1_k_ = fft_.fac1D()*(*sol->pLin(i, 0));
//            fft_.inverse( &tmp1_k_, &rho_); //  Real space rho
//            energy_u += rho_.abs().sum();
//            energy_b += divu_.abs().sum();
        }
        
        
    }
    double mpi_send_buff[num_to_mpi] = {energy_u,energy_b,AM_u,AM_b,diss_u,diss_b};
    double mpi_receive_buff[num_to_mpi];
    
    // Put the everything on processor 0
    mpi_.SumReduce_doub(mpi_send_buff,mpi_receive_buff,num_to_mpi);

    if (tv->reynolds_save_Q()){
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
    }
    
    // Currently, TimeVariables is set to save on root process, may want to generalize this at some point (SumReduce_doub is also)
    if (mpi_.my_n_v() == 0) {
        ////////////////////////////////////////////
        ///// All this is only on processor 0  /////
        double divfac=L(2)/totalN2_;
        
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
        double divfavMF = L(2)/( NZfull()*NZfull() );
        
        energy_MB = ( (*(sol->pMF(0))).abs2().sum() + (*(sol->pMF(1))).abs2().sum() )*divfavMF;
        energy_MU = ( (*(sol->pMF(2))).abs2().sum() + (*(sol->pMF(3))).abs2().sum() )*divfavMF;
        
        AM_MB = divfavMF*( (*sol->pMF(0))*(sol->pMF(1)->conjugate()) ).real().sum();
        AM_MU = divfavMF*( (*sol->pMF(2))*(sol->pMF(3)->conjugate()) ).real().sum();
        
        diss_MB = (divfavMF*eta_)*( -(K->kz2*sol->pMF(0)->abs2()).sum() - (K->kz2*sol->pMF(1)->abs2()).sum() );
        diss_MU = (divfavMF*nu_)*( -(K->kz2*sol->pMF(2)->abs2()).sum() - (K->kz2*sol->pMF(3)->abs2()).sum() );
        
        
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
        
        // Set etaK_times_kmax in TimeVariables - This is then used in output for monitoring spatial resolution as simulation progresses
        // Use kz_max for this scale (sometimes have high x resolution)
        // Calculate on zero and broadcast
        tv->etaK_times_kmax[0] = pow(nu_*nu_*nu_/(diss_MU+diss_u),0.25)*(2*PI/L(2)*(N(2)/3));
        tv->etaK_times_kmax[1] = pow(eta_*eta_*eta_/(diss_MB+diss_b),0.25)*(2*PI/L(2)*(N(2)/3));
        
        
        if (tv->reynolds_save_Q()) {
            //////////////////////////////////////
            //   Reynolds stress
            //  Saves the stress on each of Bx, By, Ux, Uy in this order. Each is normalized by the rms mean-field amplitude
            double* rey_point = tv->current_reynolds();
            
            if (!save_full_reynolds_){
                // Measuring Channel modes - keep B and U stresses
                // By driving/damping
                rey_point[0] = (Bx_drive_*sol->pMF(0)->conjugate()).real().sum()/sqrt(sol->pMF(0)->abs2().sum());
                // By driving/damping
                rey_point[1] = (By_drive_*sol->pMF(1)->conjugate()).real().sum()/sqrt(sol->pMF(1)->abs2().sum());
                // U Reynolds and Maxwell stresses
                // Ux driving/damping
                rey_point[2] = (Ux_drive_*sol->pMF(2)->conjugate()).real().sum()/sqrt(sol->pMF(2)->abs2().sum());
                // Uy driving/damping
                rey_point[3] = (Uy_drive_*sol->pMF(3)->conjugate()).real().sum()/sqrt(sol->pMF(3)->abs2().sum());
            } else {
                // SAVE FULL NONLINEAR STRESS DATA
                fft_.inverse( &Bx_drive_, &reynolds_z_tmp_);
                for (int jj=0; jj<NZfull(); ++jj) {
                    rey_point[jj] = real(reynolds_z_tmp_(jj))*fft_.fac1D();
                }
                fft_.inverse( &By_drive_, &reynolds_z_tmp_);
                for (int jj=0; jj<NZfull(); ++jj) {
                    rey_point[jj+NZfull()] = real(reynolds_z_tmp_(jj))*fft_.fac1D();
                }
                fft_.inverse( &Ux_drive_, &reynolds_z_tmp_);
                for (int jj=0; jj<NZfull(); ++jj) {
                    rey_point[jj+2*NZfull()] = real(reynolds_z_tmp_(jj))*fft_.fac1D();
                }
                fft_.inverse( &Uy_drive_, &reynolds_z_tmp_);
                for (int jj=0; jj<NZfull(); ++jj) {
                    rey_point[jj+3*NZfull()] = real(reynolds_z_tmp_(jj))*fft_.fac1D();
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
inline void CompMHDIso_FixedMF::assign_laplacians_(int i, double t, bool need_inverse){
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



void CompMHDIso_FixedMF::print_noise_range_(){
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
double CompMHDIso_FixedMF::Calculate_CFL(const solution *sol)  {
    // Returns CFL/dt to calculate dt - in this case CFL/dt = kmax By + q + kmax*cs
    
    // Mean fields have been previously calculated, so may as well use them
    double Bxmax = sqrt(MBx_.abs2().maxCoeff());
    double Bymax = sqrt(MBy_.abs2().maxCoeff());
    double Uxmax = sqrt(MUx_.abs2().maxCoeff());
    double Uymax = sqrt(MUy_.abs2().maxCoeff());
    // CFL
    return kmax*(Bymax + Bxmax + Uxmax + Uymax) + q_ + kmax*sqrt(cs20_);
    
}



//////////////////////////////////////////////////////////
//   DRIVING NOISE
void CompMHDIso_FixedMF::DrivingNoise(double t, double dt, solution *sol) {
    // NZ() is now dealiased NZ!!!
    // So - drive all of NZ and Nyquist frequency is not included

    
    // Adds noise onto linear part of solution
    for (int i=0; i<Dimxy(); ++i) {
        
        if (K->ky_index[i] != 0){ // ky=0 dealt with seperately. kx=ky=0 missed automatically
            
            ///////////////////////////
            // Noise multipliers
            assign_laplacians_(i,t,1); // Assign laplacian's 
            
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
            
            // In the compressible case, since we want to clean the divergence from the noise before adding it into u, generate noise into tmp1_k_, tmp2_k_ etc., then clean divergence, then add into u
            // ky != 0 modes are completely random
            // Ux
//            tmp1_k_.setConstant(0.0);
//            tmp2_k_.setConstant(0.0);
//            tmp3_k_.setConstant(0.0);
//            dcmplx* ePoint = tmp1_k_.data();
//            for (int jj=0; jj<NZ(); ++jj) {
//                if (drive_cond_pnt[jj] && !drive_only_magnetic_fluctuations_)
//                    ePoint[jj] = noise_multfac*dcmplx(ndist_(mt_), ndist_(mt_));
//            }
//            // Uy
//            ePoint = tmp2_k_.data();
//            for (int jj=0; jj<NZ(); ++jj) {
//                if (drive_cond_pnt[jj] && !drive_only_magnetic_fluctuations_)
//                    ePoint[jj] = noise_multfac*dcmplx(ndist_(mt_), ndist_(mt_));
//            }
//            // Uz
//            ePoint = tmp3_k_.data();
//            for (int jj=0; jj<NZ(); ++jj) {
//                if (drive_cond_pnt[jj] && !drive_only_magnetic_fluctuations_)
//                    ePoint[jj] = noise_multfac*dcmplx(ndist_(mt_), ndist_(mt_));
//            }
//            
//            // Clean divergence from noise and add onto U
//            tmp4_k_ = ilapFtmp_*(kxctmp_*tmp1_k_ + kyctmp_*tmp2_k_ + K->kz*tmp3_k_);
//            *(sol->pLin(i,1) ) += tmp1_k_ - kxctmp_*tmp4_k_;
//            *(sol->pLin(i,2) ) += tmp2_k_ - kyctmp_*tmp4_k_;
//            *(sol->pLin(i,3) ) += tmp3_k_ - K->kz*tmp4_k_;
            
            
            // TEMPORARY VERSION - actually, maybe this is better anyway?
            // This version calculates vy and vz noise from zeta, so as to compare with incompressible version
            // ky != 0 modes are completely random
            // U
            tmp1_k_.setZero();
            tmp2_k_.setZero();
            tmp3_k_.setZero();
            tmp4_k_.setZero();
            dcmplx* ePoint = tmp1_k_.data();
            for (int jj=0; jj<NZ(); ++jj) {
                if (drive_cond_pnt[jj] && !drive_only_magnetic_fluctuations_)
                    ePoint[jj] = multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
            }
            // zeta
            ePoint = tmp2_k_.data();
            for (int jj=0; jj<NZ(); ++jj) {
                if (drive_cond_pnt[jj] && !drive_only_magnetic_fluctuations_)
                    ePoint[jj] = multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
            }

            *(sol->pLin(i,1) ) += tmp1_k_;
            *(sol->pLin(i,2) ) += ( (-kyctmp_*kxctmp_)*tmp1_k_ +  K->kz*tmp2_k_ )*ilap2tmp_;
            *(sol->pLin(i,3) ) += ( -kxctmp_*K->kz*tmp1_k_ -  kyctmp_*tmp2_k_ )*ilap2tmp_;
            
            // Same for the magnetic field
            // b
            ePoint = tmp3_k_.data();
            for (int jj=0; jj<NZ(); ++jj) {
                if (drive_cond_pnt[jj] && !drive_only_velocity_fluctuations_)
                    ePoint[jj] += multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
            }
            // eta
            ePoint = tmp4_k_.data();
            for (int jj=0; jj<NZ(); ++jj) {
                if (drive_cond_pnt[jj] && !drive_only_velocity_fluctuations_)
                    ePoint[jj] += multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
            }
            *(sol->pLin(i,4) ) += tmp3_k_;
            *(sol->pLin(i,5) ) += ( (-kyctmp_*kxctmp_)*tmp3_k_ +  K->kz*tmp4_k_ )*ilap2tmp_;
            *(sol->pLin(i,6) ) += ( -kxctmp_*K->kz*tmp3_k_ -  kyctmp_*tmp4_k_ )*ilap2tmp_;
        }
        
        // Don't drive kx=ky=0 modes (mean fields) or ky=kz=0 modes (non-invertible).These are both taken care of in laplacian
        // MODIFIED -- NOW WILL DRIVE ky=kz=0 MODES
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
            
            assign_laplacians_(i,t,1); // Assign laplacian's - no inverses
            
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
            // Use tmp1_k_ etc. to store noise so that divergence constraint can be enforced on the noise. I.e., tmp1_k and tmp2_k act just like u and zeta until right at the end, after all of the message passing
            tmp1_k_.setZero();
            tmp2_k_.setZero();
            tmp3_k_.setZero();
            tmp4_k_.setZero();
            dcmplx* noise_buff_dat_ = noise_buff_.data();
            dcmplx* ePoint = tmp1_k_.data();
            for (int jj=1; jj<NZ(); ++jj){
                if (drive_cond_pnt[jj] && !drive_only_magnetic_fluctuations_){
                    noise_buff_dat_[jj-1]=multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                    ePoint[jj] += noise_buff_dat_[jj-1];
                }
            }
            // zeta
            noise_buff_dat_ += noise_buff_len_;
            ePoint = tmp2_k_.data();
            for (int jj=1; jj<NZ(); ++jj){
                if (drive_cond_pnt[jj] && !drive_only_magnetic_fluctuations_){
                    noise_buff_dat_[jj-1]=multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                    ePoint[jj] += noise_buff_dat_[jj-1];
                }
            }
            // Assign u and zeta noise to ux, uy, uz
            *(sol->pLin(i,1) ) += tmp1_k_;
            *(sol->pLin(i,2) ) += ( (-kyctmp_*kxctmp_)*tmp1_k_ +  K->kz*tmp2_k_ )*ilap2tmp_;
            *(sol->pLin(i,3) ) += ( -kxctmp_*K->kz*tmp1_k_ -  kyctmp_*tmp2_k_ )*ilap2tmp_;
            
            
            // b
            noise_buff_dat_ += noise_buff_len_;
            ePoint = tmp3_k_.data();
            for (int jj=1; jj<NZ(); ++jj){
                if (drive_cond_pnt[jj] && !drive_only_velocity_fluctuations_){
                    noise_buff_dat_[jj-1]=multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                    ePoint[jj] += noise_buff_dat_[jj-1];
                }
            }
            // eta
            noise_buff_dat_ += noise_buff_len_;
            ePoint = tmp4_k_.data();
            for (int jj=1; jj<NZ(); ++jj){
                if (drive_cond_pnt[jj] && !drive_only_velocity_fluctuations_){
                    noise_buff_dat_[jj-1]=multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                    ePoint[jj] += noise_buff_dat_[jj-1];
                }
            }
            // Assign b and eta noise to bx, by, bz
            *(sol->pLin(i,4) ) += tmp3_k_;
            *(sol->pLin(i,5) ) += ( (-kyctmp_*kxctmp_)*tmp3_k_ +  K->kz*tmp4_k_ )*ilap2tmp_;
            *(sol->pLin(i,6) ) += ( -kxctmp_*K->kz*tmp3_k_ -  kyctmp_*tmp4_k_ )*ilap2tmp_;
            
            // Send data
            mpi_.Send_dcmplx(noise_buff_.data(), num_noise_vars_*noise_buff_len_, *(K->i_sp), *(K->i_tosend));
            
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
            // Assign k's again, since need for conversion to ux etc.
            assign_laplacians_(i,t,1);
            
            // Receive data from matching call
            mpi_.Recv_dcmplx(noise_buff_.data(), num_noise_vars_*noise_buff_len_, *(K->i_fp), *(K->i_from));
            
            // Flip noise around
            for (int nV=0; nV<num_noise_vars_; ++nV) {
                noise_buff_.segment(nV*noise_buff_len_, noise_buff_len_).reverseInPlace();
            }
            

            tmp1_k_.setZero();
            tmp2_k_.setZero();
            tmp3_k_.setZero();
            tmp4_k_.setZero();
            { // Loop through each variable
                int nV=0;
                tmp1_k_.segment(1,NZ()-1) += noise_buff_.segment(nV*noise_buff_len_, noise_buff_len_).conjugate();
                ++nV;
                tmp2_k_.segment(1,NZ()-1) += noise_buff_.segment(nV*noise_buff_len_, noise_buff_len_).conjugate();
                ++nV;
                tmp3_k_.segment(1,NZ()-1) += noise_buff_.segment(nV*noise_buff_len_, noise_buff_len_).conjugate();
                ++nV;
                tmp4_k_.segment(1,NZ()-1) += noise_buff_.segment(nV*noise_buff_len_, noise_buff_len_).conjugate();
            }
            // Assign u and zeta noise to ux, uy, uz
            *(sol->pLin(i,1) ) += tmp1_k_;
            *(sol->pLin(i,2) ) += ( (-kyctmp_*kxctmp_)*tmp1_k_ +  K->kz*tmp2_k_ )*ilap2tmp_;
            *(sol->pLin(i,3) ) += ( -kxctmp_*K->kz*tmp1_k_ -  kyctmp_*tmp2_k_ )*ilap2tmp_;
            // Assign b and eta noise to bx, by, bz
            *(sol->pLin(i,4) ) += tmp3_k_;
            *(sol->pLin(i,5) ) += ( (-kyctmp_*kxctmp_)*tmp3_k_ +  K->kz*tmp4_k_ )*ilap2tmp_;
            *(sol->pLin(i,6) ) += ( -kxctmp_*K->kz*tmp3_k_ -  kyctmp_*tmp4_k_ )*ilap2tmp_;
            
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
            assign_laplacians_(i,t,1); // Assign laplacian's - no inverses
            
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
            
            noise_buff_.setZero();
            // Make some noise - add to current k value
            tmp1_k_.setZero();
            tmp2_k_.setZero();
            tmp3_k_.setZero();
            tmp4_k_.setZero();
            dcmplx* noise_buff_dat_ = noise_buff_.data();
            dcmplx* ePoint = tmp1_k_.data();
            for (int jj=1; jj<NZ(); ++jj){
                if (drive_cond_pnt[jj] && !drive_only_magnetic_fluctuations_){
                    noise_buff_dat_[jj-1]=multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                    ePoint[jj] += noise_buff_dat_[jj-1];
                }
            }
            // zeta
            noise_buff_dat_ += noise_buff_len_;
            ePoint = tmp2_k_.data();
            for (int jj=1; jj<NZ(); ++jj){
                if (drive_cond_pnt[jj] && !drive_only_magnetic_fluctuations_){
                    noise_buff_dat_[jj-1]=multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                    ePoint[jj] += noise_buff_dat_[jj-1];
                }
            }
            // Assign u and zeta noise to ux, uy, uz
            *(sol->pLin(i,1) ) += tmp1_k_;
            *(sol->pLin(i,2) ) += ( (-kyctmp_*kxctmp_)*tmp1_k_ +  K->kz*tmp2_k_ )*ilap2tmp_;
            *(sol->pLin(i,3) ) += ( -kxctmp_*K->kz*tmp1_k_ -  kyctmp_*tmp2_k_ )*ilap2tmp_;
            
            
            // b
            noise_buff_dat_ += noise_buff_len_;
            ePoint = tmp3_k_.data();
            for (int jj=1; jj<NZ(); ++jj){
                if (drive_cond_pnt[jj] && !drive_only_velocity_fluctuations_){
                    noise_buff_dat_[jj-1]=multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                    ePoint[jj] += noise_buff_dat_[jj-1];
                }
            }
            // eta
            noise_buff_dat_ += noise_buff_len_;
            ePoint = tmp4_k_.data();
            for (int jj=1; jj<NZ(); ++jj){
                if (drive_cond_pnt[jj] && !drive_only_velocity_fluctuations_){
                    noise_buff_dat_[jj-1]=multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                    ePoint[jj] += noise_buff_dat_[jj-1];
                }
            }
            // Assign b and eta noise to bx, by, bz
            *(sol->pLin(i,4) ) += tmp3_k_;
            *(sol->pLin(i,5) ) += ( (-kyctmp_*kxctmp_)*tmp3_k_ +  K->kz*tmp4_k_ )*ilap2tmp_;
            *(sol->pLin(i,6) ) += ( -kxctmp_*K->kz*tmp3_k_ -  kyctmp_*tmp4_k_ )*ilap2tmp_;
            
            ////////////////////////////////////////////
            // Find matching K value and use the same noise
            i = *(K->i_loc);
            
            // Assign k's again, since need for conversion to ux etc.
            assign_laplacians_(i,t,1);
            
            // Flip noise around
            for (int nV=0; nV<num_noise_vars_; ++nV) {
                noise_buff_.segment(nV*noise_buff_len_, noise_buff_len_).reverseInPlace();
            }
            
            
            tmp1_k_.setZero();
            tmp2_k_.setZero();
            tmp3_k_.setZero();
            tmp4_k_.setZero();
            { // Loop through each variable
                int nV=0;
                tmp1_k_.segment(1,NZ()-1) += noise_buff_.segment(nV*noise_buff_len_, noise_buff_len_).conjugate();
                ++nV;
                tmp2_k_.segment(1,NZ()-1) += noise_buff_.segment(nV*noise_buff_len_, noise_buff_len_).conjugate();
                ++nV;
                tmp3_k_.segment(1,NZ()-1) += noise_buff_.segment(nV*noise_buff_len_, noise_buff_len_).conjugate();
                ++nV;
                tmp4_k_.segment(1,NZ()-1) += noise_buff_.segment(nV*noise_buff_len_, noise_buff_len_).conjugate();
            }
            // Assign u and zeta noise to ux, uy, uz
            *(sol->pLin(i,1) ) += tmp1_k_;
            *(sol->pLin(i,2) ) += ( (-kyctmp_*kxctmp_)*tmp1_k_ +  K->kz*tmp2_k_ )*ilap2tmp_;
            *(sol->pLin(i,3) ) += ( -kxctmp_*K->kz*tmp1_k_ -  kyctmp_*tmp2_k_ )*ilap2tmp_;
            // Assign b and eta noise to bx, by, bz
            *(sol->pLin(i,4) ) += tmp3_k_;
            *(sol->pLin(i,5) ) += ( (-kyctmp_*kxctmp_)*tmp3_k_ +  K->kz*tmp4_k_ )*ilap2tmp_;
            *(sol->pLin(i,6) ) += ( -kxctmp_*K->kz*tmp3_k_ -  kyctmp_*tmp4_k_ )*ilap2tmp_;
            
            
            // Update iterators
            ++(K->i_tosend);
            ++(K->i_loc);
        }
        
#endif
        ////////////////////////////////////////////////////////////////////
        
    }
    
}
//////////////////////////////////////////////////////////

