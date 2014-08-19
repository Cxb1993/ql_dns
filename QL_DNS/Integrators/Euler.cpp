//
//  Euler.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/19/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "Euler.h"

Euler::Euler(double t0, double dt, Model * mod) :
dt_(dt), nxy_(mod->Dimxy()),
model_(mod)
{
    solution sol_rhs_ = new solution(mod);
}

Euler::~Euler()
{ delete sol_rhs_;
}


// Step forward in time
int Euler::Step(double t, solution * sol){
    
    model_->rhs(t, 0, sol ,sol_rhs_);
    
    for (int i=0; i<nxy_; ++i) {
        &(sol->pLin(i)) += &( sol_rhs_->pLin(i) )*dt_;
    }
}