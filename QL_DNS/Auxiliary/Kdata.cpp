//
//  Kdata.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/22/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "Kdata.h"

Kdata::Kdata(Model *mod, MPIdata *mpi){
    //   First form entire array, then take bit for each processor

    //////////////////////////////////////////
    // Define kx and ky arrays - length is Nx*Ny
    // Leave out the (zero,zero) frequncy
    dcmplx *kx_tmp, *ky_tmp; // kx and ky - each of length nxy
    ;
    kx_tmp = new dcmplx[ mod->Dimxy_full() ];
    ky_tmp = new dcmplx[ mod->Dimxy_full() ];
    ky_index_full = new int[ mod->Dimxy_full() ];
    kx_index_full = new int[ mod->Dimxy_full() ];

    kxLfac = 2*PI/( mod->L(0) );
    kyLfac = 2*PI/( mod->L(1) );
    kzLfac = 2*PI/( mod->L(2) );
    
    // Include only positive ky values!
    int ny=mod->N(1)/2, nx=mod->N(0)/2; // Nx and Ny over 2!
    for (int i=0; i<nx; ++i) {
        for (int j=0; j<ny; ++j) {
            kx_tmp[j+ny*i]=dcmplx(0,i*kxLfac);
            ky_tmp[j+ny*i]=dcmplx(0,j*kyLfac);
            kx_index_full[j+ny*i]=i;
            ky_index_full[j+ny*i]=j;
        }
    }
    // Leave out the Nyquist frequncy in kx
    for (int i=nx+1; i<2*nx; ++i) {
        for (int j=0; j<ny; ++j) {
            kx_tmp[j+ny*(i-1)]=dcmplx(0,(-2*nx+i)*kxLfac);
            ky_tmp[j+ny*(i-1)]=dcmplx(0,j*kyLfac);
            kx_index_full[j+ny*(i-1)]=-2*nx+i;
            ky_index_full[j+ny*(i-1)]=j;
        }
    }
    ////////////////////////////////////////////////
    // Take bit belonging to each processor
    int nxy = mod->Dimxy();
    kx = new dcmplx[ nxy ];
    ky = new dcmplx[nxy];
    ky_index = new int[ nxy ];
    kx_index = new int[ nxy ];
    for (int i=0; i<nxy ;  ++i){
        int k_i = i +  mpi->minxy_i(); // k index
        kx[i] = kx_tmp[k_i];
        ky[i] = ky_tmp[k_i];
        ky_index[i] = ky_index_full[k_i];
        kx_index[i] = kx_index_full[k_i];
    }

    
    ///////////////////////////////////////////
    // kz and kz^2 arrays
    // kz is length NZ eigen array
    kz = dcmplxVec( mod->NZ() );
    kz2 = doubVec( mod->NZ() );
    int nzl=(mod->NZ() - 1)/2; // Loop nz
    for (int i=0; i<nzl+1; ++i)
        kz(i) = dcmplx(0,i*kzLfac );
    for (int i=nzl+1; i<2*nzl+1; ++i)
        kz(i) = dcmplx(0,(-2*nzl+i-1)*kzLfac );
    
    kz2 = (kz*kz).real();



    
/////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////
    //   EXTRA BITS
    // lap2 arrays
    // Wasting memory a bit here, since all ky are stored on all processors. To improve, would be better to have a more general parallelization
    // Maximum ky index
    int number_of_ky = ny;
    // Quick error check
    if (mod->Dimxy_full()%number_of_ky!=0)
        std::cout << "Warning, something funny going on in Define_Lap2_Arrays_" << std::endl;
    // Assign data to arrays
    lap2 = new doubVec[ number_of_ky ];
    ilap2 = new doubVec[ number_of_ky ];
    
    
    for (int i=0; i<number_of_ky;  ++i){
        
        // Form Laplacian
        doubVec lap2tmp( mod->NZ() ),ilap2tmp( mod->NZ() );
        double kyr_tmp=ky_tmp[i].imag();
        lap2tmp = -kyr_tmp*kyr_tmp + kz2;
        ilap2tmp = 1/lap2tmp;
        
        // Fix up zero bits
        if (kyr_tmp==0 ) {
            lap2tmp(0)=0;
            // Avoid infinities
            ilap2tmp(0)=1;
        }
        
        // Assign to lap2 and ilap2
        lap2[i] = lap2tmp;
        ilap2[i] = ilap2tmp;
        
        
    }

/////////////////////////////////////////////////
    // Arrays for communications to create real noise
    int *ky0pkx = new int[nx-1]; // Locations of ky=0 in array, kx>0
    int *ky0nkx = new int[nx-1]; // Locations of ky=0 in array, kx<0
    for (int i=0; i<2*nx-2; ++i){
        if (i<nx-1)
            ky0pkx[i] = (i+1)*ny;
        if (i>=nx-1)
            ky0nkx[2*(nx-1)-i-1] = (i+1)*ny;
    }
    
    for (int i=0; i<nx-1; ++i) {
        if (ky0pkx[i]/nxy == mpi->my_n_v()) {
            // Add value to vector if it is in ky0 array
            match_kx_local.push_back(ky0pkx[i]%nxy);
            // Add where it receives this from
            match_kx_from_p.push_back(ky0nkx[i]/nxy);
            match_kx_from.push_back(ky0nkx[i]%nxy);  // Don't need this but use as a tag to make sure the comm matches up
        }
        if (ky0nkx[i]/nxy == mpi->my_n_v()){
            // Add value to _tosend vectors
            match_kx_tosend.push_back(ky0nkx[i]%nxy);
            match_kx_tosend_p.push_back(ky0pkx[i]/nxy);
        }
    }
    
    // SAMPLE CODE TO PERFORM NECESSARY SENDS AND RECEIVES
    
//    i_tosend = match_kx_tosend.begin();
//    i_loc = match_kx_local.begin();
//    i_from = match_kx_from.begin();
//    i_fp = match_kx_from_p.begin();
//    i_sp = match_kx_tosend_p.begin();
//    // Sending
//    while (i_sp < match_kx_tosend_p.end()) {
//        double tosend[2] = {kx[*i_tosend].imag(),ky[*i_tosend].imag()};
//        MPI_Send(&tosend, 2, MPI_DOUBLE, *i_sp, *i_tosend, MPI_COMM_WORLD);
//        ++ i_tosend;
//        ++ i_sp;
//    }
//    // Recieving
//    while (i_fp < match_kx_from_p.end()) {
//        double torec[2];
//        MPI_Recv(torec, 2, MPI_DOUBLE, *i_fp, *i_from, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        cout << torec[0] << " " << kx[*i_loc].imag() << ", " << torec[1] << " " << ky[*i_loc].imag() << endl;
//        
//        ++i_fp;
//        ++i_loc;
//        ++i_from;
//    }
    
    delete[] kx_tmp; // Could make these members if it ends up being required
    delete[] ky_tmp;
    delete[] ky0pkx;
    delete[] ky0nkx;
}


Kdata::~Kdata() {
    delete[] kx;
    delete[] ky;
    delete[] ky_index;
    delete[] kx_index;
    delete[] ky_index_full;
    delete[] kx_index_full;
    
    delete[] lap2;
    delete[] ilap2;
}