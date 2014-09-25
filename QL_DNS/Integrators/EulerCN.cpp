//
//  EulerCN.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/21/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "EulerCN.h"


EulerCN::EulerCN(double t0, Inputs& SP, Model * mod) :
nxy_(mod->Dimxy()), nz_(mod->NZ()),
nMF_(mod->num_MFs()),nLin_(mod->num_Lin()),
model_(mod)
{
    if (SP.dt<0)
        std::cout << "Variable time-step not supported by EulerCN integrator!"<< std::endl;
    dt_ = SP.dt; // This integrator requires a time-step to be specified in inputs!!
    
    // Solution RHS (dt U = F(U) ) as returned by integrator
    sol_rhs_ = new solution(model_);
    
    
    /////////////////////////////////////////////////////////////
    // Linear operators - couldn't think of a good way to do this, just doubVec array
    
    linearOp_fluct_ = new doubVec*[nxy_];
    linearOp_fluct_old_ = new doubVec*[nxy_];
    for (int i=0; i<nxy_; ++i) {
        linearOp_fluct_[i] = new doubVec[nLin_];
        linearOp_fluct_old_[i] = new doubVec[nLin_];
        for (int j=0; j<nLin_; ++j) {
            linearOp_fluct_[i][j] = doubVec(nz_);
            linearOp_fluct_old_[i][j] = doubVec(nz_);
        }
    }
    denom = doubVec(nz_);
    
    doubVec* linearOp_MF_ = new doubVec[nMF_];
    linearOp_MF_NLCo_ = new doubVec[nMF_];
    linearOp_MF_LinCo_ = new doubVec[nMF_];
    for (int i=0; i<nMF_; ++i) {
        linearOp_MF_[i] = doubVec(nz_);
        linearOp_MF_NLCo_[i] = doubVec(nz_);
        linearOp_MF_LinCo_[i] = doubVec(nz_);
    }
    model_->linearOPs_Init(t0, linearOp_fluct_old_, linearOp_MF_);
    
    // Pre-assign the MF linear operators
    // ASSUMES MF OPERATOR CONSTANT IN TIME!
    for (int i=0; i<nMF_; ++i) {
        linearOp_MF_NLCo_[i] = dt_/(1.0-dt_/2*linearOp_MF_[i]);
        linearOp_MF_LinCo_[i] = (1.0+dt_/2*linearOp_MF_[i])/(1.0-dt_/2*linearOp_MF_[i]);
    }
    delete[] linearOp_MF_;
}

EulerCN::~EulerCN() {
    delete sol_rhs_;
    for (int i=0; i<nxy_; ++i) {
        delete[] linearOp_fluct_[i];
        delete[] linearOp_fluct_old_[i];
    }
    delete[] linearOp_fluct_;
    delete[] linearOp_fluct_old_;
    delete[] linearOp_MF_LinCo_;
    delete[] linearOp_MF_NLCo_;
    
}


///////////////////////////////////////////
//   STEP
// Step forward in time
double EulerCN::Step(double t, solution * sol){
    
    model_->rhs(t, dt_, sol ,sol_rhs_, linearOp_fluct_);
    
    // Can't think of a nice way to isolate all the details of sol from the integrator for the more complicated semi-implicit versions - annoying!
    // Thus, this depends on the details of the solution class!
    
    // Linear fluctuations
    for (int i=0; i<nxy_; ++i) {
        for (int j=0; j<nLin_; ++j){
            denom = 1/(1-dt_/2*linearOp_fluct_[i][j]);
            // Probably faster to do inside Eigen for vectorization
            (*(sol->pLin(i,j))) = dt_*denom*(*(sol_rhs_->pLin(i,j))) + (1+dt_/2*linearOp_fluct_old_[i][j])*denom*(*(sol->pLin(i,j)));
        }
        
        
    }
    
    // Mean fields
    for (int i=0; i<nMF_; ++i) {
        (*sol->pMF(i)) = linearOp_MF_NLCo_[i]*(*sol_rhs_->pMF(i)) + linearOp_MF_LinCo_[i]*(*sol->pMF(i));
    }
    
    // Move variables around
    for (int i=0; i<nxy_; i++){
        for (int j=0; j<nLin_; ++j)
            linearOp_fluct_old_[i][j] = linearOp_fluct_[i][j];
    }
    
    return dt_;

    
    
}
