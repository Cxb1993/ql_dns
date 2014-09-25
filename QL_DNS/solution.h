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
#include "Auxiliary/fftwPlans.h"


class solution{
public:
    // Constructor - takes a model
    solution(Model* eqs);
    ~solution();
    
    ////  REQUIRED METHODS  -  pLin, pMF
    // Pointers to parts of solutions
    dcmplxVec* pLin(int i,int Vnum) const {return linear_field_[i]+Vnum;}; //  Fluctuation variables
    dcmplxVec* pMF(int i) const {return mean_field_+i;}; //  Fluctuation variables
    
    // Initial conditions
    void Initial_Conditions(fftwPlans &fft);
        
    
private:
    // Sizes - nz, nxy (dimension in x,y on each processor), nMF (number of mean fields), nLin (number of linear fields)
    const int nz_, nz_full_, nxy_, nMF_, nLin_;
    // Data storage
    dcmplxVec **linear_field_;// Fluctuations - unlike CE2 store fields individually. With container class don't gain much from joining these together
    dcmplxVec *mean_field_; // Mean fields
    
    
};


#endif /* defined(__QL_DNS__solution__) */





