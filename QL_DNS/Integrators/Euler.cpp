//
//  Euler.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/19/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "Euler.h"

Euler::Euler(double t0, double dt, Model * mod) :
dt_(dt), nxy_(mod->Dimxy()), nz_(mod->NZ()),
nMF_(mod->num_MFs()),nLin_(mod->num_Lin()),
model_(mod)
{
    // Solution RHS (dt U = F(U) ) as returned by integrator
    sol_rhs_ = new solution(model_);
    
    /////////////////////////////////////////////////////////////
    // Linear operators - couldn't think of a good way to do this, just doubVec array
    linearOp_fluct_ = new doubVec[nxy_];
    for (int i=0; i<nxy_; ++i) {
        linearOp_fluct_[i] = doubVec(nLin_*nz_);
    }
    linearOp_MF_ = new doubVec[nMF_];
    for (int i=0; i<nMF_; ++i) {
        linearOp_MF_[i] = doubVec(nz_);
    }
    model_->linearOPs_Init(t0, linearOp_fluct_, linearOp_MF_);
}

Euler::~Euler() {
    delete sol_rhs_;
    delete[] linearOp_fluct_;
    delete[] linearOp_MF_;

}


// Step forward in time
void Euler::Step(double t, solution * sol){
    
    model_->rhs(t, 0, sol ,sol_rhs_,linearOp_fluct_);
        
    // Can't think of a nice way to isolate all the details of sol from the integrator for the more complicated semi-implicit versions - annoying!
    
    // Crappy and unessessary but can't be bothered changing
    for (int i=0; i<nxy_; ++i) {
        dLin_ = (*(sol->pLin(i))).data();
        dLin_rhs_ = (*(sol_rhs_->pLin(i))).data();
        lin_op_ = linearOp_fluct_[i].data();
        
        for (int j=0; j<nz_*nLin_; ++j) {
            dLin_[j] += dt_*dLin_rhs_[j] + dt_*lin_op_[j]*dLin_[j];
        }
    }
    
    // Mean fields
    for (int i=0; i<nMF_; ++i) {
        (*sol->pMF(i)) += dt_*(*sol_rhs_->pMF(i)) + dt_*linearOp_MF_[i]*(*sol->pMF(i));
    }
}