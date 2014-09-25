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
    void calculatePlans( int NZnd, int nz  );
    
    // Size
    double fac1D() const {return ftfac_;};// Factor to divide by in inverse
    
    /////////////////////////////////////////////////////
    // Functions for exectuting the various FFTs, given pointers to data
    // Including dealias or padding with zeros
    //  DOES NOT INCLUDE MULTIPLICATION BY 1/NZ IN INVERSE!
    void forward(dcmplxVec *invec, dcmplxVec *outvec){
#ifndef EIGEN_NO_DEBUG // Check lengths are correct
        if (invec->size() != NZnd_ || outvec->size() != nz_) {
            std::cout << "Sizes of vectors in forward transform are incompatible!!" << std::endl;
        }
#endif
        for_1D(invec->data());
        dealias(invec->data(), outvec->data());
    }
    
    void inverse(dcmplxVec *invec, dcmplxVec *outvec){
#ifndef EIGEN_NO_DEBUG // Check lengths are correct
        if (invec->size() != nz_ || outvec->size() != NZnd_) {
            std::cout << "Sizes of vectors in inverse transform are incompatible!!" << std::endl;
        }
#endif
        pad_with_zeros(invec->data(), outvec->data());
        back_1D(outvec->data());
    }
    ////////////////////////////////////////////////////
    // Basic 1-D transforms 1-D MF transforms
    void for_1D(dcmplx* MFin){
        fftw_execute_dft(t1D_for_,CAST_T0_FFTW(MFin),CAST_T0_FFTW(MFin));
    };
    void back_1D(dcmplx* MFin){
        fftw_execute_dft(t1D_back_,CAST_T0_FFTW(MFin),CAST_T0_FFTW(MFin));
    };
    /////////////////////////////////////////////////////

    
    /////////////////////////////////////////////////////
    // DEALIASING - keep public for occasional use
    // Function to pad with zeros (inverse dealais) - backward transform
    void pad_with_zeros(dcmplx* invec, dcmplx* outvec){
        // invec - length nz_  ( NZ() )
        // outvec - length NZnd_ ( N(2) )
        int i=0,j=deal_low_;
        while (i<deal_low_) {
            outvec[i] = invec[i];
            ++i;
        }
        while (j<deal_high_){
            outvec[j] = 0.0;
            ++j;
        }
        while (j<NZnd_) {
            outvec[j] = invec[i];
            ++i;++j;
        }
    };
    // For forwards transform
    void dealias(dcmplx* invec, dcmplx* outvec){
        // invec - length NZnd_  ( N(2) )
        // outvec - length nz_ ( NZ() )
        
        int i=0,j=deal_high_;
        while (i<deal_low_) {
            outvec[i] = invec[i];
            ++i;
        }
        while (j<NZnd_) {
            outvec[i] = invec[j];
            ++i;++j;
        }
    };
    
private:

    // Sizes
    int NZnd_; // Non-dealiased length
    int nz_; // Dealiased length
    double ftfac_; // 1/NZ for normalization
    int deal_low_, deal_high_; // Bounds for copying data to pad with zeros

    
    // Transforms of 1-D data - same for MF and fluctuations
    fftw_plan t1D_for_, t1D_back_;  // Forward 1-D transform
    
    bool plans_calculated_; // 1 if calculatePlans has been run
    
    
};


#endif /* defined(__QL_DNS__fftwPlans__) */
