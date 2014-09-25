//
//  Euler.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/19/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "Euler.h"

Euler::Euler(double t0, Inputs& SP, Model * mod) :
nxy_(mod->Dimxy()), nz_(mod->NZ()),
nMF_(mod->num_MFs()),nLin_(mod->num_Lin()),
model_(mod)
{
    if (SP.dt<0)
        std::cout << "Variable time-step not supported by Euler integrator!" << std::endl;
    dt_ = SP.dt; // This integrator requires a time-step to be specified in input!!
    
    // Solution RHS (dt U = F(U) ) as returned by integrator
    sol_rhs_ = new solution(model_);
    
    /////////////////////////////////////////////////////////////
    // Linear operators - couldn't think of a good way to do this, just doubVec array
    linearOp_fluct_ = new doubVec*[nxy_];
    for (int i=0; i<nxy_; ++i) {
        linearOp_fluct_[i] = new doubVec[nLin_];
        for (int Vn=0; Vn<nLin_; ++Vn) {
            linearOp_fluct_[i][Vn] = doubVec(nz_);
        }
    }
    linearOp_MF_ = new doubVec[nMF_];
    for (int i=0; i<nMF_; ++i) {
        linearOp_MF_[i] = doubVec(nz_);
    }
    model_->linearOPs_Init(t0, linearOp_fluct_, linearOp_MF_);
}

Euler::~Euler() {
    delete sol_rhs_;
    for (int i=0; i<nxy_; ++i) {
        delete[] linearOp_fluct_[i];
    }
    delete[] linearOp_fluct_;
    delete[] linearOp_MF_;

}


// Step forward in time
double Euler::Step(double t, solution * sol){
    
    model_->rhs(t, 0, sol ,sol_rhs_,linearOp_fluct_);
        
    // Can't think of a nice way to isolate all the details of sol from the integrator for the more complicated semi-implicit versions - annoying!
    
    // Crappy and unessessary but can't be bothered changing
    for (int i=0; i<nxy_; ++i) {
        for (int j=0; j<nLin_; ++j){
            // Probably faster to do inside Eigen for vectorization
            (*(sol->pLin(i,j))) = dt_*(*(sol_rhs_->pLin(i,j))) + linearOp_fluct_[i][j]*(*(sol->pLin(i,j)));
        }

    }
    
    // Mean fields
    for (int i=0; i<nMF_; ++i) {
        (*sol->pMF(i)) += dt_*(*sol_rhs_->pMF(i)) + dt_*linearOp_MF_[i]*(*sol->pMF(i));
    }
    
    return dt_;
}