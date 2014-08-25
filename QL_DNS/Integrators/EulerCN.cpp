//
//  EulerCN.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/21/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "EulerCN.h"


EulerCN::EulerCN(double t0, double dt, Model * mod) :
dt_(dt), nxy_(mod->Dimxy()), nz_(mod->NZ()),
nMF_(mod->num_MFs()),nLin_(mod->num_Lin()),
model_(mod)
{
    // Solution RHS (dt U = F(U) ) as returned by integrator
    sol_rhs_ = new solution(model_);
    
    
    /////////////////////////////////////////////////////////////
    // Linear operators - couldn't think of a good way to do this, just doubVec array
    linearOp_fluct_ = new doubVec[nxy_];
    linearOp_fluct_old_ = new doubVec[nxy_];
    for (int i=0; i<nxy_; ++i) {
        linearOp_fluct_[i] = doubVec(nLin_*nz_);
        linearOp_fluct_old_[i] = doubVec(nLin_*nz_);
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
    delete[] linearOp_fluct_;
    delete[] linearOp_fluct_old_;
    delete[] linearOp_MF_LinCo_;
    delete[] linearOp_MF_NLCo_;
    
}


///////////////////////////////////////////
//   STEP
// Step forward in time
void EulerCN::Step(double t, solution * sol){
    
    model_->rhs(t, 0, sol ,sol_rhs_, linearOp_fluct_);
    
    // Can't think of a nice way to isolate all the details of sol from the integrator for the more complicated semi-implicit versions - annoying!
    // Thus, this depends on the details of the solution class!
    
    // Linear fluctuations
    for (int i=0; i<nxy_; ++i) {
        denom = 1/(1-dt_/2*linearOp_fluct_[i]);
        // Probably faster to do inside Eigen for vectorization
        (*(sol->pLin(i))) = dt_*denom*(*(sol_rhs_->pLin(i))) + (1+dt_/2*linearOp_fluct_old_[i])*denom*(*(sol->pLin(i)));
    }
    
    // Mean fields
    for (int i=0; i<nMF_; ++i) {
        (*sol->pMF(i)) = linearOp_MF_NLCo_[i]*(*sol_rhs_->pMF(i)) + linearOp_MF_LinCo_[i]*(*sol->pMF(i));
    }
    
    // Move variables around
    for (int i=0; i<nxy_; i++)
        linearOp_fluct_old_[i] = linearOp_fluct_[i];
    
    
}
