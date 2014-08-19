//
//  Euler.h
//  QL_DNS
//
//  Created by Jonathan Squire on 8/19/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __QL_DNS__Euler__
#define __QL_DNS__Euler__

#include "Integrator.h"

class Euler : public Integrator {
public:
    Euler(double t0, double dt, Model * mod);
    ~Euler();
    // Step
    int Step(double t, dcmplxVec * sol);
    
private:
    const double dt_;
    const int nxy_;
    
    // Model
    Model * model_;
    
    // Space for integrator evaluation
    solution * sol_rhs_; // Output of model
};



#endif /* defined(__QL_DNS__Euler__) */
