//
//  solution.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/10/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "solution.h"


solution::solution(Model* eqs) :
nMF_(eqs->num_MFs()), nLin_(eqs->num_Lin()), nz_(eqs->NZ()), nxy_(eqs->Dimxy())
{
    // Constructor for the solution class, contains data required for holding the QL solution

    linear_field_ = new dcmplxVec[nxy_];
    for (int i=0; i<nxy_; ++i) {
        linear_field_[i] = dcmplxVec(nLin_*nz_);
    }
    
    mean_field_ = new dcmplxVec[nMF_];
    for (int i=0; i<nMF_; ++i) {
        mean_field_[i] = dcmplxVec(nz_);
    }
    
}

solution::~solution() {
    delete[] linear_field_;
    delete[] mean_field_;
}


// Initial conditions
void solution::Initial_Conditions() {
    for (int i=0; i<nxy_; ++i) {
        linear_field_[i].setConstant(1.0);
    }
    for (int i=0; i<nMF_; ++i) {
        mean_field_[i].setConstant(1.0);
    }
}





