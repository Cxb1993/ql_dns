//
//  main.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/10/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "General_Definitions.h"
// Models
#include "Models/ConstantDamping.h"
// Integrators
#include "Integrators/Euler.h"
#include "Integrators/EulerCN.h"
// Auxiliary
#include "Auxiliary/MPIdata.h"
#include "Auxiliary/Input_parameters.h"


int main(int argc, char ** argv)
{
    MPIdata mpi; // Storage of MPI data
#ifdef USE_MPI_FLAG
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, mpi.total_n_p());
    MPI_Comm_rank(MPI_COMM_WORLD, mpi.my_n_p());
    MPI_Comm_size(MPI_COMM_WORLD, mpi.comm_size_p());
#endif
    
    
    ////////////////////////////////////////////////////////////////
    /////////                                           ////////////
    /////////                INPUTS                     ////////////
    std::string input_file_name;
    if (argc == 2) {
        input_file_name = argv[1]; // If an argument is passed, this specifies input file
    } else {
        input_file_name = "null";
    }
    Inputs SP( mpi , input_file_name);
    /////////                                           ////////////
    ////////////////////////////////////////////////////////////////

    double t=SP.t_start ; // Initial time

    // Create a Model
    fftwPlans fft;
    Model* fluidEqs = new ConstantDamping(SP, mpi, fft);
    
    // Initialize solution
    solution *sol = new solution( fluidEqs );
    
    // Initial condtions
    sol->Initial_Conditions();

    
    // Integrator
    Integrator *integrator = new EulerCN(t, SP.dt, fluidEqs);
    
    // Main loop
    for (int i = SP.i_start+1; i < SP.nsteps+1; i++) {
        
        std::stringstream out;
        out << (*(sol->pLin(3))).transpose() << std::endl << (*(sol->pMF(0))).transpose() << std::endl;
        mpi.print1(out.str());

        
        integrator->Step(t, sol);
        t += SP.dt;
        
    }
    
    
    delete fluidEqs;
    delete sol;
    delete integrator;
#ifdef USE_MPI_FLAG
    MPI_Finalize();
#endif
    return 0;
}

