//
//  MHD_FullQlin.h
//  QL_DNS
//
//  Created by Jonathan Squire on 9/16/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __QL_DNS__MHD_FullQlin__
#define __QL_DNS__MHD_FullQlin__


#include "Model.h"
#include "../solution.h"
// Auxiliary classes
#include "../Auxiliary/Kdata.h"
#include "../Auxiliary/fftwPlans.h"
#include "../Auxiliary/Input_parameters.h"
#include "../Auxiliary/MPIdata.h"


class MHD_FullQlin : public Model {
public:
    MHD_FullQlin(const Inputs& sp, MPIdata& mpi, fftwPlans& fft);
    ~MHD_FullQlin();
    
    // Equations
    void rhs(const double t, const double dt_lin,
              solution * SolIn, solution * SolOut, doubVec **linOpFluct);
    
    // Linear operator initialization
    void linearOPs_Init(double t0, doubVec **linOpFluct, doubVec *linOpMF);
    
    
    // Sizes
    int Dimxy() const {return mpi_.nxy();}; // Size of C array on given MPI process
    int num_MFs() const {return numMF_;};  // Number of mean fields
    int num_Lin() const {return numLin_;}; // Number of fluctuating fields
    
    int num_reynolds_saves(){return num_reynolds_saves_;};
    
    //////////////////////////////////////////////////////////
    //   DRIVING NOISE
    void DrivingNoise(double t, double dt, solution *sol);
    // Helpful for generating initial conditions
    void ChangeNoiseRange(double fnoise, double kmin, double kmax){
        f_noise_=fnoise;
        ndist_ = boost::random::normal_distribution<double>(0,f_noise_/sqrt(2));
        noise_range_[0] = kmin*kmin;
        noise_range_[1] = kmax*kmax;};
    //////////////////////////////////////////////////////////
    
    //////////////////////////////////////////////////////////
    //   Remapping
    void ShearingBox_Remap(double qt, solution *sol);
    
    //////////////////////////////////////////////////////////////////
    //////  AUXILIARY FUNCTIONS OPERATING ON SOLUTION   //////////////
    void Calc_Energy_AM_Diss(TimeVariables* tv, double t, const solution* sol);
    //////////////////////////////////////////////////////////////////
    
    //////////////////////////
    // CFL number
    double Calculate_CFL(const solution *sol);
    double kmax;
    
private:
    const std::string equations_name; // Name of model
    // Number of fields for this model
    const int numMF_, numLin_;
    // Physical parameters - don't set as const
    double nu_, eta_; // dissipation
    const int viscosity_order_; // Order of the Laplacian operator for dissipation
    double f_noise_; // Driving noise
    const double q_; // q
    // Mean magnetic field
    const double B0z_;
    
    bool QL_YN_; // Flag for QL or not
    
    bool dont_drive_ky0_modes_;
    bool drive_only_velocity_fluctuations_;
    bool drive_only_magnetic_fluctuations_;
    bool save_full_reynolds_;
    bool turn_off_fluct_Lorentz_force_;
    int L_mult_;
    int num_reynolds_saves_;
    
    const double Omega_; // Rotation level, 0 for pure shear, 2/3*q for Keplerian
    
    // fft
    fftwPlans &fft_;
    
    // mpi
    MPIdata &mpi_;
    
    // Useful array sizes
    long totalN2_;
    double mult_noise_fac_;
    
    // Assign laplacians - just for convenience
    void assign_laplacians_(int i,double t,bool need_inverse);
    
    // Random number generation
    boost::random::mt19937 mt_;
    boost::random::normal_distribution<double> ndist_;
    dcmplxVec noise_buff_;
    int noise_buff_len_; // Length of noise (no need to include dealiased bits)
    // High/low noise cutoff
    double noise_range_[2];
    Eigen::Matrix<bool,Eigen::Dynamic,1> drive_condition_; // Store driving in z as you loop  though x,y
    void print_noise_range_();
    
    
    // Reynolds stresses
    dcmplxVec Bx_drive_,By_drive_, Ux_drive_, Uy_drive_;
    dcmplxVec reynolds_z_tmp_;// Temporary, size NZfull()
    dcmplxVec reynolds_stress_MPI_send_, reynolds_stress_MPI_receive_;
    int num_to_mpi_send_; // Number of terms to send for Reynolds stress per MF variable
    double reynolds_stress_ft_fac; // 1/(nx*ny)^2
    
    ////////////////////////////////////////////////////
    //               TEMPORARY VARIABLES              //
    doubVec lapFtmp_, lap2tmp_; // Laplacians - nice to still have lap2
    doubVec ilapFtmp_, ilap2tmp_; // Inverse Laplacians
    
    // Temps for evaulation of main equations
    // u and b and derivatives - Real space
    dcmplxVec *u_, *b_; // (uy,uz,dx{ux,uy,uz},dy{ux,uy,uz}) (and b)
    dcmplxVec fluct_y_,fluct_z_; // In Fourier space
    dcmplxVec tmp1_k_,tmp2_k_,tmp3_k_; // Temps before transform (size NZ() )
    dcmplxVec tmp1_z_;
    
    dcmplx kxctmp_, kyctmp_;
    double kxtmp_,kytmp_;
    int ind_ky_;
    
    dcmplxVec MBy_,dzMBy_,MBx_,dzMBx_; // Mean fields (real space)
    dcmplxVec MUy_,dzMUy_,MUx_,dzMUx_; // Mean fields (real space)

    
};

#endif /* defined(__QL_DNS__MHD_FullQlin__) */
