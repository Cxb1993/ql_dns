//
//  fftwPlans.h
//  QL_DNS
//
//  Created by Jonathan Squire on 8/23/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __QL_DNS__fftwPlans__
#define __QL_DNS__fftwPlans__

#include "../General_Definitions.h"


// FFTW plan storage (and execution) class
class fftwPlans {
    // Stores 1-D and 2-D transforms for the mean fields and Ckl
public:
    // Constructor
    fftwPlans(): plans_calculated_(0) {};
    // Destructor
    ~fftwPlans();
    // No Copy constructor, but should only have one instance anyway
    // Setup plans
    void calculatePlans( int NZ );
    
    /////////////////////////////////////////////////////
    // Functions for exectuting the various FFTs, given pointers to data
    // 1-D MF transforms
    void for_1D(dcmplx* MFin){
        fftw_execute_dft(t1D_for_,CAST_T0_FFTW(MFin),CAST_T0_FFTW(MFin));
    };
    void back_1D(dcmplx* MFin){
        fftw_execute_dft(t1D_back_,CAST_T0_FFTW(MFin),CAST_T0_FFTW(MFin));
    };
    
    /////////////////////////////////////////////////////
    
    
private:

    // Transforms of 1-D data - same for MF and fluctuations
    fftw_plan t1D_for_, t1D_back_;  // Forward 1-D transform
    
    bool plans_calculated_; // 1 if calculatePlans has been run
    
};


#endif /* defined(__QL_DNS__fftwPlans__) */
