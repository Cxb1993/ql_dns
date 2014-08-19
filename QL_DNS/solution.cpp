//
//  solution.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/10/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "solution.h"

solution::solution(Model* eqs) {
    // Constructor for the solution class, contains data required for holding the QL solution
    nz_ = eqs->NZ();
    dimxy_ = eqs->Dimxy_full();
    
    linear_field_ = new dcmplxVec[dimxy_];
    for (int i=0; i<dimxy_; ++i) {
        linear_field_[i] = dcmplxVec(nz_);
    }
    
}

solution::~solution() {
    delete[] linear_field_;
}


// Initial conditions
void solution::Initial_Conditions() {
    for (int i=0; i<dimxy_; ++i) {
        linear_field_[i].setConstant(1.0);
    }
}