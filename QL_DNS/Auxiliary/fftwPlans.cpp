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
void fftwPlans::calculatePlans( int NZ ) {

    ////////////////////////////////////////////////
    //      CREATE TEMPORARY EIGEN OBJECTS FOR PLANS
    //  I'm relatively sure it doesn't create fftw problems if you later delete these...
    
    // Vector temps
    dcmplxVec MF_tmp1( NZ );
    dcmplx * MF_p1 = MF_tmp1.data();
    
    
    ////////////////////////////////////////////////
    //   1D Plans
    //  Need both a 1-D plan and a 1-D*NZ plan
    //
    // Only in place transforms have been useful, could delete others
    
    
    // True 1-D plans
    t1D_for_ = fftw_plan_dft_1d(NZ, CAST_T0_FFTW(MF_p1), CAST_T0_FFTW(MF_p1),
                                FFTW_FORWARD, MY_FFTWPLAN);
    t1D_back_ = fftw_plan_dft_1d(NZ, CAST_T0_FFTW(MF_p1), CAST_T0_FFTW(MF_p1),
                                 FFTW_BACKWARD, MY_FFTWPLAN);
    
    
    // Flag to use in destructor
    plans_calculated_ = 1;
 
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
