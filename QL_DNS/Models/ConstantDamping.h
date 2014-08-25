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

// Very basic model class to use for ConstantDampinging things. Just sets solOut = solIn, so simple exponential growth
class ConstantDamping : public Model {
public:
    ConstantDamping(const Inputs& sp, MPIdata& mpi, fftwPlans& fft);
    ~ConstantDamping();
    
    // Equations
    void rhs(const double t, const double dt_lin,
             const solution * SolIn, solution * SolOut, doubVec *linOpFluct);
    
    // Linear operator initialization
    void linearOPs_Init(double t0, doubVec *linOpFluct, doubVec *linOpMF);
    
    
    // Sizes
    int Dimxy() const {return Dimxy_full();}; // Size of C array on given MPI process
    int num_MFs() const {return numMF_;};  // Number of mean fields
    int num_Lin() const {return numLin_;}; // Number of fluctuating fields

private:
    const std::string equations_name; // Name of model
    // Number of fields for this model
    const int numMF_, numLin_;
    // Physical parameters - don't set as const
    double nu_, eta_; // dissipation
    double f_noise_; // Driving noise
    
    double dampFac;
    
    // fft
    fftwPlans &fft_;
    
    // mpi
    MPIdata &mpi_;
};

#endif /* defined(__QL_DNS__ConstantDamping__) */
