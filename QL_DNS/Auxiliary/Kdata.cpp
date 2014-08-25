//
//  Kdata.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/22/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "Kdata.h"

Kdata::Kdata(Model *mod){
    //////////////////////////////////////////
    // Define kx and ky arrays - length is Nx*Ny
    // Leave out the (zero,zero) frequncy
    kx = new dcmplx[ mod->Dimxy_full() ];
    ky = new dcmplx[ mod->Dimxy_full() ];
    ky_index = new int[ mod->Dimxy_full() ];
    // Include only positive ky values!
    int ny=mod->N(1)/2, nx=mod->N(0)/2; // Nx and Ny over 2!
    for (int i=0; i<nx; ++i) {
        for (int j=0; j<ny; ++j) {
            kx[j+ny*i]=dcmplx(0,i*2*PI/( mod->L(0) ));
            ky[j+ny*i]=dcmplx(0,j*2*PI/( mod->L(1) ));
            ky_index[j+ny*i]=j;
        }
    }
    // Leave out the Nyquist frequncy in kx
    for (int i=nx+1; i<2*nx; ++i) {
        for (int j=0; j<ny; ++j) {
            kx[j+ny*(i-1)]=dcmplx(0,(-2*nx+i)*2*PI/( mod->L(0) ));
            ky[j+ny*(i-1)]=dcmplx(0,j*2*PI/( mod->L(2) ));
            ky_index[j+ny*(i-1)]=j;
        }
    }
    ////////////////////////////////////////////////
    
    ///////////////////////////////////////////
    // kz and kz^2 arrays
    // kz is length NZ eigen array
    kz = dcmplxVec( mod->N(2) );
    kz2 = doubVec( mod->N(2) );
    int nz=mod->N(2)/2;
    for (int i=0; i<nz; ++i)
        kz(i) = dcmplx(0,i*2*PI/(mod->L(2)) );
    for (int i=nz; i<2*nz; ++i)
        kz(i) = dcmplx(0,(-2*nz+i)*2*PI/(mod->L(2)) );
    
    kz2 = (kz*kz).real();

    
}


Kdata::~Kdata() {
    delete[] kx;
    delete[] ky;
    delete ky_index;
}