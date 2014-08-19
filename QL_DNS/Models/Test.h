//
//  Test.h
//  QL_DNS
//
//  Created by Jonathan Squire on 8/19/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __QL_DNS__Test__
#define __QL_DNS__Test__

#include "Model.h"

// Very basic model class to use for testing things. Just sets solOut = solIn, so simple exponential growth
class Test : public Model {
public:
    Test(int nz, int * nxy, double growthFac);
    ~Test();
    
    // Equations
    void rhs(const double t, const double dt_lin,
             const dcmplxVec * SolIn, dcmplxVec * SolOut);
    
    // Sizes
    int Dimxy() const {return Dimxy_full();}; // Size of C array on given MPI process
    int num_MFs() const {return numMF_;};  // Number of mean fields

private:
    const int numMF_;
    const double growthFac_;
};

#endif /* defined(__QL_DNS__Test__) */
