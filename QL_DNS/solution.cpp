//
//  solution.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/10/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "solution.h"


solution::solution(Model* eqs) :
nMF_(eqs->num_MFs()), nLin_(eqs->num_Lin()), nz_(eqs->NZ()), nz_full_(eqs->NZfull()), nxy_(eqs->Dimxy())
{
    // Constructor for the solution class, contains data required for holding the QL solution
    // ONLY STORES NON-DEALIASED VALUES OF SOLUTION IN FOURIER SPACE (eqs->NZ() returns this)
    // IN kz STORED AS: (0 1 2 3 ... (NZ()-1)/2 -(NZ()-1)/2  -(NZ()-1)/2+1 ... -2 -1)
    
    linear_field_ = new dcmplxVec *[nxy_];
    for (int i=0; i<nxy_; ++i) {
        linear_field_[i] = new dcmplxVec[nLin_];
        for (int j=0; j<nLin_; ++j) {
            linear_field_[i][j] = dcmplxVec(nz_);
        }
    }
    
    mean_field_ = new dcmplxVec[nMF_];
    for (int i=0; i<nMF_; ++i) {
        mean_field_[i] = dcmplxVec(nz_);
    }
    
}

solution::~solution() {
    for (int i=0; i<nxy_; ++i) {
        delete[] linear_field_[i];
    }
    delete[] linear_field_;
    delete[] mean_field_;
}


// Initial conditions
void solution::Initial_Conditions(fftwPlans &fft) {
    for (int i=0; i<nxy_; ++i) {
        for (int j=0; j<nLin_; ++j) {
            linear_field_[i][j].setConstant(0.0);
        }
        
    }
    // MEAN FIELDS
    dcmplxVec meanf_r = dcmplxVec::Constant(nz_full_, 0.1); // Real space version
    mean_field_[0].setConstant(0.0);
    fft.forward(&meanf_r, mean_field_+1);

}





