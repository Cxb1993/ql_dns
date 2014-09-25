//
//  fftwPlans.cc
//  QL_DNS
//
//  Created by Jonathan Squire on 8/23/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "fftwPlans.h"

/////////////////////////////////////////////////////
//                                                 //
//                 fftw_Plans                      //
//                                                 //
/////////////////////////////////////////////////////

// fftw plans constructor - creates the plans
void fftwPlans::calculatePlans( int NZnd, int nz )
{

    // FFTW container class for calculating and storing plans
    // Also contains delaising functions, since these are really padding with zeros before taking ffts
    
    ////////////////////////////////////////////////
    //      CREATE TEMPORARY EIGEN OBJECTS FOR PLANS
    //  I'm relatively sure it doesn't create fftw problems if you later delete these...
    
    
    NZnd_ = NZnd; // Non-dealiased length
    ftfac_ = 1.0/NZnd_;
    nz_ = nz; // Dealiased length
    // Dealiasing bounds
    deal_low_ = (nz_-1)/2+1;
    deal_high_ = NZnd_ - (nz_-1)/2;
    
    // Assign memory to calculate the plan
    dcmplx *buffer = new dcmplx[NZnd_];
    
    
    ////////////////////////////////////////////////
    //   1D Plans
    //  Need both a 1-D plan and a 1-D*NZ plan
    //
    // Only in place transforms have been useful, could delete others
    
    
    // True 1-D plans
    t1D_for_ = fftw_plan_dft_1d(NZnd_, CAST_T0_FFTW(buffer), CAST_T0_FFTW(buffer),
                                FFTW_FORWARD, MY_FFTWPLAN);
    t1D_back_ = fftw_plan_dft_1d(NZnd_, CAST_T0_FFTW(buffer), CAST_T0_FFTW(buffer),
                                 FFTW_BACKWARD, MY_FFTWPLAN);
    
    
    // Flag to use in destructor
    plans_calculated_ = 1;
    
    delete[] buffer;
    // NB: Might be slightly faster to use advanced interface for all the Eigen data at once, but this would require changing a lot of stuff..
}

// Destructor
fftwPlans::~fftwPlans() {
    ///////////////////////////////////////////////////////////////////////
    // This deletes data that has been created in fftwPlans::calculatePlans
    
    // std::cout << "Destroying " << std::endl;
    if (plans_calculated_) {
        // 1-D plans
        fftw_destroy_plan(t1D_for_);
        fftw_destroy_plan(t1D_back_);
    }
    
    
    
    fftw_cleanup();
    
}
//                                                 //
/////////////////////////////////////////////////////
