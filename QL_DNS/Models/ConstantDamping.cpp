//
//  ConstantDamping.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/19/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "ConstantDamping.h"

ConstantDamping::ConstantDamping(const Inputs& sp, MPIdata& mpi, fftwPlans& fft) :
dampFac(0.5),
equations_name("ConstantDamping"),
numMF_(1), numLin_(4),
f_noise_(sp.f_noise), nu_(sp.nu), eta_(sp.eta),
q_(sp.q),
dont_drive_ky0_modes_(0),// If true, no driving ky=0 modes
Model(sp.NZ, sp.NXY , sp.L), // Dimensions - stored in base
mpi_(mpi), // MPI data
fft_(fft) // FFT data
{
    // Setup MPI
    mpi_.Split_NXY_Grid( Dimxy_full() ); // Work out MPI splitting
    // Assign K data
    K = new Kdata(this, &mpi_); // Stores all the K data
    // Fourier transform plans
    fft_.calculatePlans( N(2),NZ() );
    
    // Random generator
    mt_ = boost::random::mt19937(1.0 + mpi_.my_n_v());  // Seed is 1.0, could change
    ndist_ = boost::random::normal_distribution<double>(0,f_noise_/2); // Normal distribution, standard deviation f_noise_ (factor 2 is since it's complex)
    noise_buff_ = dcmplxVec(num_Lin()*(NZ()-1));// NZ-1 since ky=kz=0 mode is not driven
    
    // Sizes of various arrays used for normalizing things
    totalN2_ = N(0)*N(1)*N(2); // Total number of grid points
    totalN2_ = totalN2_*totalN2_; // Squared
    mult_noise_fac_ = 1.0/(16*32*32); // Defined so that consistent with (16, 32, 32) results
    mult_noise_fac_ = mult_noise_fac_*mult_noise_fac_;
    // NZ^2 for normalizing  mean field energy
    nz2_ = N(2)*N(2);
    
    
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


}


ConstantDamping::~ConstantDamping(){
    delete K;
}

void ConstantDamping::rhs(const double t, const double dt_lin,
               const solution * SolIn, solution * SolOut,doubVec **linOpFluct) {
    
    for (int i=0; i<Dimxy(); ++i) {
        // Full Loop containing all the main work equation
        ///////////////////////////////////////
        ///// MAIN EQUATIONS
        
        assign_laplacians_(i, t, 1);
        for (int Vn=0; Vn<num_Lin(); ++Vn) {
            *( SolOut->pLin(i,Vn) ) = -(*( SolIn->pLin(i,Vn) ))*dampFac;
        }
        
        
    
        ////////////////////////////////////////
        //////   LINEAR PART
        // Need to re-evaluate laplacian, since different time.
        kxtmp_=kxtmp_ + q_*dt_lin*kytmp_;
        
        linOpFluct[i][0].setZero();
        linOpFluct[i][1].setZero();
        linOpFluct[i][2].setZero();
        linOpFluct[i][3].setZero();

//        linOpFluct[i][0] = nu_*((-kxtmp_*kxtmp_-kytmp_*kytmp_)+ K->kz2);
//        linOpFluct[i][1] = linOpFluct[i][0];
//        linOpFluct[i][2]  = (eta_/nu_)*linOpFluct[i][0];// Saves some calculation
//        linOpFluct[i][3] = linOpFluct[i][2];
        
        ////////////////////////////////////////
    
        
    }
    for (int i=0; i<num_MFs(); ++i) {
        (*( SolOut->pMF(i) )).setZero();;
    }
    

}

void ConstantDamping::linearOPs_Init(double t0, doubVec **linOpFluct, doubVec *linOpMF){
    // Fluctuations
    for (int i=0; i<Dimxy(); ++i) {
        for (int j=0; j<num_Lin(); ++j) {
            linOpFluct[i][j].setZero();
        }
    }
    // Mean fields
    for (int i=0; i<num_MFs(); ++i) {
        linOpMF[i].setZero();
    }
}


//////////////////////////////////////////////////////////
//   Remapping
// Remapping for the shearing box
// Using the Lithwick method of continuously remapping wavenumbers when they become too large
void ConstantDamping::ShearingBox_Remap(double qt, solution *sol){
    int nx = N(0)-1;
    
    double kxt;
    for (int i=0; i<Dimxy(); ++i) {
        kxt = K->kx[i].imag() + qt*K->ky[i].imag();
        // If kx has grown too large
        if (kxt > nx/2*K->kxLfac) {
            // Put kx back in correct range
            K->kx[i] = K->kx[i] - dcmplx(0,nx*K->kxLfac);
            // zero out solution
            for (int Vn=0; Vn<num_Lin(); ++Vn) {
                (sol->pLin(i, Vn))->setZero();
            }
        }
    }
}
//////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////
//    ENERGY/MOMENTUM ETC.
void ConstantDamping::Calc_Energy_AM_Diss(TimeVariables* tv, double t, const solution *sol) {
    // Energy, angular momentum and dissipation of the solution MFin and Cin
    // TimeVariables class stores info and data about energy, AM etc.
    // t is time
    
    // OUTPUT: energy[1] and [2] contain U and B mean field energies (energy[1]=0 for this)
    // energy[3] and [4] contain u and b fluctuating energies.
    
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
        
        if (tv->energy_save_Q()){
            //////////////////////////////////////
            //     ENERGY
            // Use Qkl_tmp_ for Mkl to save memory
            lap2tmp_ = lapFtmp_*ilap2tmp_; // lap2tmp_ just a convenient storage
            
            // (NB: ilap2tmp_ is negative but want abs, just use -ilap2)
            energy_u += mult_fac*(  lap2tmp_*(( *(sol->pLin(i,0)) ).abs2()) - ilap2tmp_*(( *(sol->pLin(i,1)) ).abs2())  ).sum();
            energy_b += mult_fac*(  lap2tmp_*(( *(sol->pLin(i,2)) ).abs2()) - ilap2tmp_*(( *(sol->pLin(i,3)) ).abs2())  ).sum();
            //////////////////////////////////////
        }
        if (tv->AngMom_save_Q()){
            //////////////////////////////////////
            //   Angular momentum
            //
            //////////////////////////////////////
        }
        if (tv->dissip_save_Q()){
            //////////////////////////////////////
            //     DISSIPATION
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
        
        energy_MB = ( (*(sol->pMF(0))).abs2().sum() )/nz2_;
       
        
        AM_MB = 0;
        AM_MU = 0;
        
        diss_MB = 0;
        diss_MU = 0;
        
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
    
    // Save the data
    tv->Save_Data(t);
    
}


//////////////////////////
// CFL number
double ConstantDamping::Calculate_CFL() const {
    // Returns CFL/dt to calculate dt - in this case don't know
    
    // Pointlessly simple
    return dampFac*4;
    
}



//////////////////////////////////////////////////////////
//   AUXILIARY FUNCTIONS
inline void ConstantDamping::assign_laplacians_(int i, double t, bool need_inverse){
    ind_ky_ = K->ky_index[i];
    // Form Laplacians using time-dependent kx
    kyctmp_ = K->ky[i];
    kytmp_ = kyctmp_.imag(); // It is useful to have both complex and double versions
    kxctmp_ = K->kx[i] + q_*t*kyctmp_;
    kxtmp_ = kxctmp_.imag();
    
    lap2tmp_ = K->lap2[ind_ky_];
    lapFtmp_ = -kxtmp_*kxtmp_ + lap2tmp_;
    if (need_inverse){
        ilapFtmp_ = 1/lapFtmp_;
        ilap2tmp_ = K->ilap2[ind_ky_];
    }
}
//////////////////////////////////////////////////////////








//////////////////////////////////////////////////////////
//   DRIVING NOISE
void ConstantDamping::DrivingNoise(double t, double dt, solution *sol) {
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
        int nz = NZ()-1; // For noise
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
            noise_buff_dat_ = noise_buff_.data()+nz;
            ePoint = (*(sol->pLin(i,1) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                noise_buff_dat_[jj-1]=multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                ePoint[jj] += noise_buff_dat_[jj-1];
            }
            // b
            noise_buff_dat_ = noise_buff_.data()+2*nz;
            ePoint = (*(sol->pLin(i,2) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                noise_buff_dat_[jj-1]=multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                ePoint[jj] += noise_buff_dat_[jj-1];
            }
            // eta
            noise_buff_dat_ = noise_buff_.data()+3*nz;
            ePoint = (*(sol->pLin(i,3) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                noise_buff_dat_[jj-1]=multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                ePoint[jj] += noise_buff_dat_[jj-1];
            }
            
            // Send data
            mpi_.Send_dcmplx(noise_buff_.data(), num_Lin()*nz, *(K->i_sp), *(K->i_tosend));
            
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
            mpi_.Recv_dcmplx(noise_buff_.data(), num_Lin()*nz, *(K->i_fp), *(K->i_from));
            
            // Flip noise around
            for (int nV=0; nV<num_Lin(); ++nV) {
                noise_buff_.segment(nV*nz, nz).reverseInPlace();
            }
            
            //        cout << "(kx,ky)=(" << K->kx[i] <<"," << K->ky[i] << ")" << endl;
            //        cout << noise_buff_.segment(0, nz).transpose().conjugate() <<endl<<endl;
            
            // Add to solution
            for (int nV=0; nV<num_Lin(); ++nV) {
                (*(sol->pLin(i,nV) )).segment(1,nz) += noise_buff_.segment(nV*nz, nz).conjugate();
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
            lap2tmp_(NZ()/2)=0.0; // Or Nyquist
            
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
            noise_buff_dat_ = noise_buff_.data()+nz;
            ePoint = (*(sol->pLin(i,1) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                noise_buff_dat_[jj-1]=multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                ePoint[jj] += noise_buff_dat_[jj-1];
            }
            // b
            noise_buff_dat_ = noise_buff_.data()+2*nz;
            ePoint = (*(sol->pLin(i,2) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                noise_buff_dat_[jj-1]=multU_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                ePoint[jj] += noise_buff_dat_[jj-1];
            }
            // eta
            noise_buff_dat_ = noise_buff_.data()+3*nz;
            ePoint = (*(sol->pLin(i,3) )).data();
            for (int jj=1; jj<NZ(); ++jj){
                noise_buff_dat_[jj-1]=multZeta_pnt[jj]*dcmplx(ndist_(mt_), ndist_(mt_));
                ePoint[jj] += noise_buff_dat_[jj-1];
            }
            
            
            // Find matching K value and use the same noise
            i = *(K->i_loc);
            
            // Flip noise around
            for (int nV=0; nV<num_Lin(); ++nV) {
                noise_buff_.segment(nV*nz, nz).reverseInPlace();
            }
            
            // Add to solution
            for (int nV=0; nV<num_Lin(); ++nV) {
                (*(sol->pLin(i,nV) )).segment(1,nz) += noise_buff_.segment(nV*nz, nz).conjugate();
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

