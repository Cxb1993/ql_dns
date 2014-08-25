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

// Class to store K data
class Kdata {
public:
    Kdata(Model *mod);
    ~Kdata();
    
    // While this is apparently bad design practice, it seems most convenient here to not use public function methods
    dcmplx *kx, *ky; // kx and ky - each of length nxy
    int *ky_index; // Index for ky (just ky/(2*pi/Ly) )
    dcmplxVec kz; // kz - more convenient as eigen object
    doubVec kz2; // kz^2, double for speed
    
};

#endif /* defined(__QL_DNS__Kdata__) */
