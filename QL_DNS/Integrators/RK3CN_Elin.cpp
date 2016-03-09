//
//  RK3CN_Elin_Elin.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 9/17/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

// Third order RK time integrator as in Lundblahd et al. 1992
// Checked error(dt) using energy and it seems to be perfectly 3rd order in time
// This version uses exp(q0*dt) for the linear part, which seems more stable for higher order hyper-resistivities
#include "RK3CN_Elin.h"


RK3CN_Elin::RK3CN_Elin(double t0, Inputs& SP, Model * mod) :
nxy_(mod->Dimxy()), nz_(mod->NZ()),
nMF_(mod->num_MFs()),nLin_(mod->num_Lin()),
model_(mod),
a0_(8.0/15.0),a1_(5.0/12.0),a2_(3.0/4.0),b1_(-17.0/60.0),b2_(-5.0/12.0),
p1_(8.0/15.0),p2_(2.0/3.0),// Integrator parameters
step_count_(0), dt_mean_(0.0)
{
    if (SP.dt < 0 ){
        variable_dt_ = 1;
        CFLnum_ = SP.CFL;
        dtmax_ = -SP.dt;
    } else {
        dt_ = SP.dt;
        variable_dt_ = 0;
    }
    
    
    // Solution RHS (dt U = F(U) ) as returned by integrator
    sol_rhs_ = new solution(model_);
    sol_rhs2_ = new solution(model_);
    
    
    /////////////////////////////////////////////////////////////
    // Linear operators - couldn't think of a good way to do this, just doubVec array
    
    // This version doesn't bother too much with the time dependence of the linear operators, just uses the final version in the exponential
    linearOp_fluct_ = new doubVec*[nxy_];
    for (int i=0; i<nxy_; ++i) {
        linearOp_fluct_[i] = new doubVec[nLin_];
        for (int j=0; j<nLin_; ++j) {
            linearOp_fluct_[i][j] = doubVec(nz_);
        }
    }
    exponL = doubVec(nz_);
    
    linearOp_MF_ = new doubVec[nMF_]; // Can't preassign since dt variable
    for (int i=0; i<nMF_; ++i) {
        linearOp_MF_[i] = doubVec(nz_);
    }
    // Don't think I actually need this step anymore, but keeps things consistent
    model_->linearOPs_Init(t0, linearOp_fluct_, linearOp_MF_);
    
}

RK3CN_Elin::~RK3CN_Elin() {
    delete sol_rhs_;
    delete sol_rhs2_;
    for (int i=0; i<nxy_; ++i) {
        delete[] linearOp_fluct_[i];
    }
    delete[] linearOp_fluct_;
    delete[] linearOp_MF_;
}


///////////////////////////////////////////
//   STEP
// Step forward in time
double RK3CN_Elin::Step(double t, solution * sol){
    
    if (variable_dt_) {
        dt_ = CFLnum_/model_->Calculate_CFL(sol);
        if (dt_ > dtmax_)
            dt_ = dtmax_;
        
    }
    dt_mean_ += dt_;
    ++step_count_;
    
    
    // Still not SURE linear part is correct on this, should check...
    //////////////////////////////////////////////////////
    //   STEP 1
    linC_ =dt_*(a0_);
    NLC1_ = a0_*dt_;
    lin_dt_ = dt_*p1_;
    
    model_->rhs(t, lin_dt_, sol ,sol_rhs_, linearOp_fluct_);
    
    // Can't think of a nice way to isolate all the details of sol from the integrator for the more complicated semi-implicit versions - annoying!
    // Thus, this depends on the details of the solution class!
    
    // Linear fluctuations
    for (int i=0; i<nxy_; ++i) {// Loop over kx ky
        
        for (int j=0; j<nLin_; ++j){ // Loop over variables
            exponL = (linC_*linearOp_fluct_[i][j]).exp();
            // Probably faster to do inside Eigen for vectorization
            (*(sol->pLin(i,j))) = exponL*( NLC1_*(*(sol_rhs_->pLin(i,j))) + (*sol->pLin(i,j)) );
        }
    }
    
    // Mean fields
    for (int i=0; i<nMF_; ++i) {
        exponL = (linC_*linearOp_MF_[i]).exp();
        (*sol->pMF(i)) = exponL*( NLC1_*(*sol_rhs_->pMF(i)) + (*sol->pMF(i)) );
    }
    
    
    //    END - STEP 1
    //////////////////////////////////////////////////////
    
    
    
    //////////////////////////////////////////////////////
    //   STEP 2 - Updates sol
    linC_ =dt_*(a1_+b1_);
    NLC1_ = a1_*dt_;
    NLC2_ = b1_*dt_;
    lin_dt_ = dt_*(p2_-p1_);
    
    model_->rhs(t+dt_*p1_, lin_dt_, sol ,sol_rhs2_, linearOp_fluct_);
    

    // Linear fluctuations
    for (int i=0; i<nxy_; ++i) {// Loop over kx ky
        
        for (int j=0; j<nLin_; ++j){ // Loop over variables
            exponL = (linC_*linearOp_fluct_[i][j]).exp();
            // Probably faster to do inside Eigen for vectorization
            (*(sol->pLin(i,j))) = exponL*( NLC1_*(*(sol_rhs2_->pLin(i,j))) + NLC2_*(*(sol_rhs_->pLin(i,j))) + (*sol->pLin(i,j)) );
        }
    }
    
    // Mean fields
    for (int i=0; i<nMF_; ++i) {
        exponL = (linC_*linearOp_MF_[i]).exp();
        (*sol->pMF(i)) = exponL*( NLC1_*(*sol_rhs2_->pMF(i)) + NLC2_*(*sol_rhs_->pMF(i)) + (*sol->pMF(i)) );
    }


    //    END - STEP 2
    //////////////////////////////////////////////////////
    
    
    
    
    //////////////////////////////////////////////////////
    //   STEP 3 - Updates sol
    linC_ =dt_*(a2_+b2_);
    NLC1_ = a2_*dt_;
    NLC2_ = b2_*dt_;
    lin_dt_ = dt_*(1.0-p2_);
    
    
    model_->rhs(t + dt_*p2_, lin_dt_, sol ,sol_rhs_, linearOp_fluct_);
    
    // Linear fluctuations
    for (int i=0; i<nxy_; ++i) {// Loop over kx ky
        
        for (int j=0; j<nLin_; ++j){ // Loop over variables
            exponL = (linC_*linearOp_fluct_[i][j]).exp();
            // Probably faster to do inside Eigen for vectorization
            (*(sol->pLin(i,j))) = exponL*( NLC1_*(*(sol_rhs_->pLin(i,j))) + NLC2_*(*(sol_rhs2_->pLin(i,j))) + (*sol->pLin(i,j)) );
        }
    }
    
    // Mean fields
    for (int i=0; i<nMF_; ++i) {
        exponL = (linC_*linearOp_MF_[i]).exp();
        (*sol->pMF(i)) = exponL*( NLC1_*(*sol_rhs_->pMF(i)) + NLC2_*(*sol_rhs2_->pMF(i)) + (*sol->pMF(i)) );
    }
    

    
    //    END - STEP 3
    //////////////////////////////////////////////////////
    
    return dt_;
    
    
    
}