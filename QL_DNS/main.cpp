//
//  main.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/10/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "General_Definitions.h"
//#include "solution.h"
// Models
#include "Models/Test.h"
// Integrators
#include "Integrators/Euler.h"


int main(int argc, const char * argv[])
{
    double dt = 0.01;
    double tfinal = 10.0;
    double t = 0;
    int nz = 16;
    int nxy[2] = {2,2};
    
    int nsteps = tfinal/dt;
    // Create a Model
    Model* fluidEqs = new Test(nz,nxy,0.1);
    
//    // Initialize solution
//    solution sol( fluidEqs );
//    
//    // Initial condtions
//    sol.Initial_Conditions();
//    
//    // Integrator
//    Integrator *integrator;// = new Euler(t, dt, fluidEqs);
//    
//    // Main loop
//    for (int i=0; i<nsteps; ++i) {
//        integrator->Step(t, &sol);
//        t += dt;
//    }
//    
    
    delete fluidEqs;
//    delete integrator;
    return 0;
}

