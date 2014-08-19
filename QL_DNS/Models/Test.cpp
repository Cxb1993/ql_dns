//
//  Test.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/19/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "Test.h"

Test::Test(int nz, int * nxy, double growthFac) :
// Constructor for simple exponential growth for testing of class structure
numMF_(0),
Model(nz,nxy),
growthFac_(growthFac)
{ 
}


Test::~Test(){
}

void Test::rhs(const double t, const double dt_lin,
               const dcmplxVec * SolIn, dcmplxVec * SolOut) {
    // Exponential growth
    for (int i=0; i<Dimxy(); ++i) {
//        &( SolOut->pLin(i) ) = &( SolIn->pLin(i) )*growthFac_;
        SolOut[i] = growthFac_*SolIn[i];
    }
}
