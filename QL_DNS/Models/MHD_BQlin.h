//
//  MHD_BQlin.h
//  QL_DNS
//
//  Created by Jonathan Squire on 9/16/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __QL_DNS__MHD_BQlin__
#define __QL_DNS__MHD_BQlin__


#include "Model.h"
#include "../solution.h"
// Auxiliary classes
#include "../Auxiliary/Kdata.h"
#include "../Auxiliary/fftwPlans.h"
#include "../Auxiliary/Input_parameters.h"
#include "../Auxiliary/MPIdata.h"

// Basic MHD Quasi-linear class
// Includes interaction of fluctuations with By and that's all
class MHD_BQlin : public Model {
public:
    MHD_BQlin(const Inputs& sp, MPIdata& mpi, fftwPlans& fft);
    ~MHD_BQlin();
    
    // Equations
    void rhs(const double t, const double dt_lin,
             const solution * SolIn, solution * SolOut, doubVec **linOpFluct);
    
    // Linear operator initialization
    void linearOPs_Init(double t0, doubVec **linOpFluct, doubVec *linOpMF);
    
    
    // Sizes
    int Dimxy() const {return mpi_.nxy();}; // Size of C array on given MPI process
    int num_MFs() const {return numMF_;};  // Number of mean fields
    int num_Lin() const {return numLin_;}; // Number of fluctuating fields
    
    
    //////////////////////////////////////////////////////////
    //   DRIVING NOISE
    void DrivingNoise(double t, double dt, solution *sol);
    //////////////////////////////////////////////////////////
    
    //////////////////////////////////////////////////////////
    //   Remapping
    void ShearingBox_Remap(double qt, solution *sol);
    
    //////////////////////////////////////////////////////////////////
    //////  AUXILIARY FUNCTIONS OPERATING ON SOLUTION   //////////////
    void Calc_Energy_AM_Diss(TimeVariables& tv, double t, const solution* sol);
    //////////////////////////////////////////////////////////////////
    
    //////////////////////////
    // CFL number
    double Calculate_CFL() const;
    double kmax;
    
private:
    const std::string equations_name; // Name of model
    // Number of fields for this model
    const int numMF_, numLin_;
    // Physical parameters - don't set as const
    double nu_, eta_; // dissipation
    double f_noise_; // Driving noise
    const double q_; // q
    
    bool QL_YN_; // Flag for QL or not
    
    bool dont_drive_ky0_modes_;
    
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
    
    // Reynolds stresses
    dcmplxVec bzux_m_uzbx_c_,bzuy_m_uzby_c_;
    dcmplxVec reynolds_z_tmp_;// Temporary, size NZfull()
    dcmplxVec reynolds_stress_MPI_send_, reynolds_stress_MPI_receive_;
    int num_to_mpi_send_; // Number of terms to send for Reynolds stress per MF variable
    double reynolds_stress_ft_fac; // 1/(nx*ny)^2
    
    ////////////////////////////////////////////////////
    //               TEMPORARY VARIABLES              //
    doubVec lapFtmp_, lap2tmp_; // Laplacians - nice to still have lap2
    doubVec ilapFtmp_, ilap2tmp_; // Inverse Laplacians
    
    // Temps for evaulation of main equations
    dcmplxVec uy_,uz_,by_,bz_;
    dcmplxVec tmp1_z_,tmp2_z_,tmp3_z_; // Temps after transform (size NZfull() )
    dcmplxVec tmp1_k_,tmp2_k_,tmp3_k_; // Temps before transform (size NZ() )
    
    dcmplx kxctmp_, kyctmp_;
    double kxtmp_,kytmp_;
    int ind_ky_;
    
    dcmplxVec By_,dzBy_,dzdzBy_; // Mean fields (real space)
    
};

#endif /* defined(__QL_DNS__MHD_BQlin__) */
