//
//  solution.h
//  QL_DNS
//
//  Created by Jonathan Squire on 8/10/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __QL_DNS__solution__
#define __QL_DNS__solution__

// General container class for the solution
// Mainly just to make it easier to pass everything around

// Stores data for fluctuations, mean-fields and the linear operators

#include "General_Definitions.h"
#include "Models/Model.h"

class solution{
public:
    // Constructor - takes a model
    solution(Model* eqs);
    ~solution();
    
    // Pointers to parts of solutions
    dcmplxVec* pLin(int i){return linear_field_+i;};
    
    // Initial condition
    void Initial_Conditions();
    
private:
    // Sizes
    int nz_, dimxy_;
    // Data storage
    dcmplxVec* linear_field_;
    
};

#endif /* defined(__QL_DNS__solution__) */
