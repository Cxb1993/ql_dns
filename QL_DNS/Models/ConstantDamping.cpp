//
//  ConstantDamping.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/19/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "ConstantDamping.h"

ConstantDamping::ConstantDamping(const Inputs& sp, MPIdata& mpi, fftwPlans& fft) :
dampFac(1.0),
equations_name("ConstantDamping"),
numMF_(1), numLin_(1),
f_noise_(sp.f_noise), nu_(sp.nu), eta_(sp.eta),
Model(sp.NZ, sp.NXY , sp.L), // Dimensions - stored in base
mpi_(mpi), // MPI data
fft_(fft) // FFT data
{
    
    // Assign K data 
    K = new Kdata(this); // Stores all the K data
    // Setup MPI
    mpi_.Split_NXY_Grid( Dimxy_full() ); // Work out MPI splitting
    // Fourier transform plans
    fft_.calculatePlans( NZ() );
    
}


ConstantDamping::~ConstantDamping(){
    delete K;
}

void ConstantDamping::rhs(const double t, const double dt_lin,
               const solution * SolIn, solution * SolOut,doubVec *linOpFluct) {
    // Exponential growth
    for (int i=0; i<Dimxy(); ++i) {
        *( SolOut->pLin(i) ) = *( SolIn->pLin(i) )*dampFac/2;
    }
    for (int i=0; i<num_MFs(); ++i) {
        *( SolOut->pMF(i) ) = *( SolIn->pMF(i) )*dampFac/2;
    }
    
    // linear part is zero
    for (int i=0; i<Dimxy(); ++i) {
        linOpFluct[i].setConstant(dampFac/2);
    }
}

void ConstantDamping::linearOPs_Init(double t0, doubVec *linOpFluct, doubVec *linOpMF){
    for (int i=0; i<Dimxy(); ++i) {
        linOpFluct[i].setConstant(dampFac/2);
    }
    for (int i=0; i<num_MFs(); ++i) {
        linOpMF[i].setConstant(dampFac/2);
    }
}