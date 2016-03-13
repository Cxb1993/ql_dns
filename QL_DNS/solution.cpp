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
    
    // Mean fields stored as (Bx, By, Ux, Uy) 
    mean_field_ = new dcmplxVec[nMF_];
    for (int i=0; i<nMF_; ++i) {
        mean_field_[i] = dcmplxVec::Zero(nz_);
    }
    mean_field_real_ = dcmplxVec(nz_full_);
    
}

solution::~solution() {
    for (int i=0; i<nxy_; ++i) {
        delete[] linear_field_[i];
    }
    delete[] linear_field_;
    delete[] mean_field_;
}




// Initial conditions
// Only used for simple cases (e.g., zero in the fluctuations). For more complex cases, save in hdf5 in Matlab and set start_from_saved_Q to 1
void solution::Initial_Conditions(Inputs &SP,fftwPlans &fft, Model *eqs, MPIdata *mpi) {
    /////////////////////////////////////////////
    //  LINEAR FIELDS
    for (int i=0; i<nxy_; ++i) {
        for (int j=0; j<nLin_; ++j) {
            linear_field_[i][j].setConstant(0.0);         
        }
    }
//    if (mpi->my_n_v()==0) linear_field_[0][1](1)=0.1;
//    linear_field_[0][0](nz_)=0.1;
//      Random Initial conditions in the linear fields - use the DrivingNoise routines
    double init_k_range[2] = {0, 100.};
    double init_amp = 0.01; // Standard deviation of the initial conditions
    // Reset model to produce driving noise like this
    eqs->ChangeNoiseRange(1.0, init_k_range[0], init_k_range[1]);
    // Add on noise to solution (which is zero)
    eqs->DrivingNoise(0.0, 1.0, this);
    double multfac = init_amp; // This has some factor in it, but figure out later based on energy. Usually want init_amp ~ 100
    if ( std::isfinite(multfac) ){
        for (int i=0; i<nxy_; ++i) {
            for (int j=0; j<nLin_; ++j) {
                linear_field_[i][j] *= multfac;
            }
        }
    } else {
        mpi->print1("Warning: Found NaN or Inf in initial conditions! Starting from zero instead");
    }
    // Reset the k range
    eqs->ChangeNoiseRange(SP.f_noise, SP.noise_range_low, SP.noise_range_high);
    //
    //////////////////////////////////////////
    
    
    ////////////////////////////////////
    //   MEAN FIELDS
    // Assign to mean fields
    
    //////////////////////////
    //  Square field/triangular velocity profiles for compressibility
//    double d4square = 3.; // d should depend on pressure to B and V ratio
//    if (eqs->num_Lin() == 4){
//        d4square = 0.01; // Revert to sine wave if incompressible
//    }
//    std::stringstream prnt;
//    prnt << "Warning: Initializing with non-sinusoidal profile with d = " << d4square<<"\n";
//    mpi->print1(prnt.str());
//    double dPonMF = 0.;
//    std::stringstream prnt;
//    prnt << "Reducing Ux and By by a factor ~ " << dPonMF<<" to account for DP0 \n";
//    mpi->print1(prnt.str());
    
    
    dcmplxVec *meanf_r = new dcmplxVec[nMF_]; // Real space version
    for (int i=0; i<nMF_; ++i) {
        meanf_r[i] = dcmplxVec(nz_full_);
    }
    doubVec zg = doubVec::LinSpaced( nz_full_, 0, 2*PI*(1.0-1.0/nz_full_) );
    
    // Decide what MF initial conditions based on SP.initial_By
    // If initial_By>0, only By nonzero, set to lowest kz mode in box
    int MODE = 1; // Mode for deterministic start
    double mult_fac;
//    if (SP.initial_By > 0.0) { // Lowest Kz mode in the box, amplitude from SP.initial_By
    int BxBy = 1; // Choice of Bx (0) or By (1)
    if (SP.initial_By >= 0.0) {//BxBy=1;} // If less than 0 use Bx
    
//        for (int i=0; i<nMF_; ++i) {
//            if (i==BxBy) { // Amplitude specified here, start By 10* larger than other(s)
//                mult_fac = fabs(SP.initial_By);
//                for (int k=0; k<nz_full_; ++k) {
//                    meanf_r[i](k) = (dcmplx) mult_fac*cos( MODE*zg(k) );
//                }
//            } else {
//                mult_fac = -0.0*SP.initial_By;
//                for (int k=0; k<nz_full_; ++k) {
//                    meanf_r[i](k) = (dcmplx) mult_fac*cos( MODE*zg(k) );
//                }
//            }
//            
//        }
        mult_fac= SP.initial_By;
        //tanh(d4square*cos( MODE*zg(k) ))/tanh(d4square);
        //sinh(d4square*sin( MODE*zg(k) ))/sinh(d4square);
        //*sqrt(1+dPonMF)
        for (int i=0; i<nMF_; ++i) {
            for (int k=0; k<nz_full_; ++k) {
                if (i==0) meanf_r[i](k) = (dcmplx) -mult_fac*cos( MODE*zg(k) );//Bx
                if (i==1) meanf_r[i](k) = (dcmplx) mult_fac*cos( MODE*zg(k) ); //By
                if (i==2) meanf_r[i](k) = (dcmplx) -sqrt(3./5.)*mult_fac*sin( MODE*zg(k) ); //Ux
                if (i==3) meanf_r[i](k) = (dcmplx) -sqrt(3./5.)*mult_fac*sin( MODE*zg(k) ); //Uy
            }
        }

    } else {// If initial_By<0, all MFs nonzero, random with specified amplitude
        for (int i=0; i<nMF_; ++i) {
            
            if (i==1) { // Amplitude specified here, start By 10* larger than other(s)
                mult_fac = -SP.initial_By;
            } else {
                mult_fac = -SP.initial_By;
            }
            meanf_r[i].real().setRandom();
            meanf_r[i].imag().setZero();
            meanf_r[i] = meanf_r[i]-meanf_r[i].mean();
            meanf_r[i] *= mult_fac;
        }
        
    }

    // Take the Fourier transform
    for (int i=0; i<nMF_; ++i) {
        fft.forward(meanf_r+i, mean_field_+i);
    }
    

    delete[] meanf_r;

}


void solution::Check_Solution(fftwPlans *fft, MPIdata *mpi) {
    dcmplxVec* check_field = mean_field_+1;
    // Check that mean fields are real
    fft->inverse(check_field, &mean_field_real_); // Check By
    
    if (mean_field_real_.real().abs().sum()/mean_field_real_.imag().abs().sum() < 1e10 &&  mean_field_real_.real().abs().sum()>1e-2){
        std::stringstream errormess;
        errormess << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        errormess << "ERROR: Mean field has developed an imaginary part!!" << std::endl;
        errormess << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        mpi->print1(errormess.str());
        ABORT;
    }
    
    // Stability
    if (check_field->abs().maxCoeff()>1e20 || mean_field_->abs().maxCoeff()>1e20) {
        std::stringstream errormess;
        errormess << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        errormess << "ERROR: Solution is very large, probably unstable!!" << std::endl;
        errormess << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        mpi->print1(errormess.str());
        ABORT;
    }
    
    //  THIS ISN'T WORKING FOR SOME REASON!!
    if (!std::isfinite(check_field->abs().sum()) || !std::isfinite(check_field->abs().sum())) {
        std::stringstream errormess;
        errormess << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        errormess << "ERROR: Solution contains NaN or Inf!!" << std::endl;
        errormess << "<<<<<<<<<<<<<<>>>>>>>>>>>>>>" << std::endl;
        mpi->print1(errormess.str());
        ABORT;
    }

}





