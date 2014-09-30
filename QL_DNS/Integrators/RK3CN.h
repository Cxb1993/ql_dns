//
//  RK3CN.h
//  QL_DNS
//
//  Created by Jonathan Squire on 9/17/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __QL_DNS__RK3CN__
#define __QL_DNS__RK3CN__

#include "Integrator.h"
#include "../solution.h"

class RK3CN : public Integrator {
public:
    RK3CN(double t0, Inputs& SP, Model * mod);
    ~RK3CN();
    // Step
    double Step(double t, solution * sol);
    
    // Average time step
    double mean_time_step() const {return dt_mean_/step_count_;};
private:
    // Time-step
    double dt_, lin_dt_;
    double CFLnum_, dtmax_,dt_mean_;
    int step_count_;
    bool variable_dt_;
    
    const int nxy_, nz_;
    const int nMF_, nLin_;
    
    // Integrator parameters
    const double a0_,a1_,a2_,b1_,b2_,p1_,p2_;
    double linC_, NLC1_, NLC2_; // Coefficients for each step
    
    
    // Model
    Model * model_;
    
    // Space for integrator evaluation
    solution * sol_rhs_, *sol_rhs2_ ; // Output of model
    
    // Space for linear operators
    doubVec **linearOp_fluct_,**linearOp_fluct_old_;
    doubVec *linearOp_MF_;
    
    doubVec denom; // Used for both mean and fluctuations
    
    //
};




#endif /* defined(__QL_DNS__RK3CN__) */


