//
//  Kdata.h
//  QL_DNS
//
//  Created by Jonathan Squire on 8/22/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __QL_DNS__Kdata__
#define __QL_DNS__Kdata__

#include "../General_Definitions.h"
#include "../Models/Model.h"
#include "MPIdata.h"

// Class to store K data
class Kdata {
public:
    Kdata(Model *mod, MPIdata *mpi);
    ~Kdata();
    
    // While this is apparently bad design practice, it seems most convenient here to not use public function methods
    dcmplx *kx, *ky; // kx and ky - each of length nxy
    int *ky_index, *kx_index; // Index for ky (just ky/(2*pi/Ly) )
    int *ky_index_full, *kx_index_full;
    dcmplxVec kz; // kz - more convenient as eigen object
    doubVec kz2; // kz^2, double for speed
    
    doubVec *lap2, *ilap2;
    
    double kxLfac,kyLfac,kzLfac;
    
    std::vector<int> match_kx_local, match_kx_from, match_kx_from_p, match_kx_tosend, match_kx_tosend_p; // Stores index of matching kx, ky=0 value for noise generation
    std::vector<int>::const_iterator i_loc, i_from, i_fp, i_tosend, i_sp;
    
};

#endif /* defined(__QL_DNS__Kdata__) */
