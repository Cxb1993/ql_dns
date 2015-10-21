//
//  ConstantDamping.h
//  QL_DNS
//
//  Created by Jonathan Squire on 8/19/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __QL_DNS__ConstantDamping__
#define __QL_DNS__ConstantDamping__

#include "Model.h"
#include "../solution.h"
// Auxiliary classes
#include "../Auxiliary/Kdata.h"
#include "../Auxiliary/fftwPlans.h"
#include "../Auxiliary/Input_parameters.h"
#include "../Auxiliary/MPIdata.h"

// Very basic model class to use for ConstantDampinging things.
// Uses standard linear MHD variables, but linear operator is just damping. Should reach equipartition of energy in all modes

class ConstantDamping : public Model {
public:
    ConstantDamping(const Inputs& sp, MPIdata& mpi, fftwPlans& fft);
    ~ConstantDamping();
    
    // Equations
    void rhs(const double t, const double dt_lin,
             const solution * SolIn, solution * SolOut, doubVec **linOpFluct);
    
    // Linear operator initialization
    void linearOPs_Init(double t0, doubVec **linOpFluct, doubVec *linOpMF);
    
    
    // Sizes
    int Dimxy() const {return mpi_.nxy();}; // Size of C array on given MPI process
    int num_MFs() const {return numMF_;};  // Number of mean fields
    int num_Lin() const {return numLin_;}; // Number of fluctuating fields
    
    int num_reynolds_saves(){return 5;};
    
    //////////////////////////////////////////////////////////
    //   DRIVING NOISE
    void DrivingNoise(double t, double dt, solution *sol);
    // Helpful for generating initial conditions
    void ChangeNoiseRange(double fnoise, double kmin, double kmax){};
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
    double Calculate_CFL(const solution *sol) ;
    
private:
    const std::string equations_name; // Name of model
    // Number of fields for this model
    const int numMF_, numLin_;
    // Physical parameters - don't set as const
    double nu_, eta_; // dissipation
    double f_noise_; // Driving noise
    const double q_; // q
    
    double dampFac;
    
    bool dont_drive_ky0_modes_;
    
    // fft
    fftwPlans &fft_;
    
    // mpi
    MPIdata &mpi_;
    
    // Useful array sizes
    long totalN2_;
    double mult_noise_fac_;
    long    nz2_;
    
    // Assign laplacians - just for convenience
    void assign_laplacians_(int i,double t,bool need_inverse);
    
    // Random number generation
    boost::random::mt19937 mt_;
    boost::random::normal_distribution<double> ndist_;
    dcmplxVec noise_buff_;
    
    
    ////////////////////////////////////////////////////
    //               TEMPORARY VARIABLES              //
    doubVec lapFtmp_, lap2tmp_; // Laplacians - nice to still have lap2
    doubVec ilapFtmp_, ilap2tmp_; // Inverse Laplacians
    dcmplxVec uy_, by_;
    
    dcmplx kxctmp_, kyctmp_;
    double kxtmp_,kytmp_;
    int ind_ky_;

};

#endif /* defined(__QL_DNS__ConstantDamping__) */
