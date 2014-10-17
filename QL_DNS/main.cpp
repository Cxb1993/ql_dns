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
#include "Models/MHD_BQlin.h"
#include "Models/MHD_FullQlin.h"
// Integrators
#include "Integrators/Euler.h"
#include "Integrators/EulerCN.h"
#include "Integrators/RK3CN.h"
// Auxiliary
#include "Auxiliary/MPIdata.h"
#include "Auxiliary/Input_parameters.h"
#include "Auxiliary/FullSave_Load.h"
#include "Auxiliary/InterfaceOutput.h"


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

    // Create a Model
    fftwPlans fft;
    Model* fluidEqs;
    std::stringstream printstr;
    printstr << "Using model: " << SP.equations_to_use << std::endl;;;
    if (SP.equations_to_use == "MHD_BQlin") {
        fluidEqs = new MHD_BQlin(SP, mpi, fft);
    } else if (SP.equations_to_use == "MHD_FullQlin") {
        fluidEqs = new MHD_FullQlin(SP, mpi, fft);
    } else if (SP.equations_to_use == "ConstantDamping") {
        fluidEqs = new ConstantDamping(SP, mpi, fft);
    } else {
        std::cout << "ERROR: no matching model found!" << std::endl;
        ABORT;
        fluidEqs = new MHD_BQlin(SP, mpi, fft); // Initialize to shut up compiler!
    }
    mpi.print1(printstr.str());
    
    
    // Initialize solution
    solution *sol = new solution( fluidEqs );
    
    // Initial condtions
    sol->Initial_Conditions(SP,fft);
    
    // Set up dump saving class (HDF5)
    FullSave_Load *sol_save = new FullSave_Load(SP, &mpi, fluidEqs); // NB: have to make this a pointer so it can be deleted before MPI_Finalize
    if (sol_save->start_from_saved()) // If starting from saved
        sol_save->LoadSolution_forRestart(sol, SP);
    
    
    double t=SP.t_start ; // Initial time
    double dt; // Time-step, set by integrator (possibly from input file)
    int step_since_TV = 0, step_since_dump = 0; // Counting steps since timevar/dump
    
    // Integrator
    Integrator *integrator = new RK3CN(t, SP, fluidEqs);

    // Set up time variables -- energy, engular momentum etc.
    TimeVariables *time_vars = new TimeVariables(SP, 4, fluidEqs->num_MFs(), 5, mpi.my_n_v());
    if (!sol_save->start_from_saved()){
        // Save initial conditions
        fluidEqs->Calc_Energy_AM_Diss(time_vars, t, sol);
        time_vars->Save_Mean_Fields(sol,  fft);
    }
    
    // Interface output and timing
    InterfaceOutput *interface = new InterfaceOutput(&mpi, time_vars, sol_save, integrator, &SP);
    
    while (t < SP.t_final  && !interface->StopSimulationQ()) {
        
        // Update solution and drive with noise
        dt = integrator->Step(t, sol);
        if (SP.f_noise != 0.0)
            fluidEqs->DrivingNoise(t,dt,sol);
        
        // Update times
        t += dt;
        SP.timevar_t += dt;
        SP.fullsave_t += dt;
        ++step_since_TV;++step_since_dump; // Steps since save
        
        
        if (SP.remapQ) // Remap at every step
            fluidEqs->ShearingBox_Remap(SP.q*t, sol);
        
        // Calculate energy, AM, etc.
        if (SP.timevar_t - SP.timvar_save_interval > -1e-8){
            fluidEqs->Calc_Energy_AM_Diss(time_vars, t, sol);
            time_vars->Save_Mean_Fields(sol,  fft);
            SP.timevar_t -= SP.timvar_save_interval;
            step_since_TV = 0;
            // Check solution
            sol->Check_Solution(&fft, &mpi);
        }
    
        // Full Solution save
        if (SP.fullsave_t - SP.fullsol_save_interval > -1e-8) {
            std::stringstream save_str;
            save_str << "Saving full solution at time " << t << std::endl;
            mpi.print1( save_str.str() );
            // Save solution
            sol_save->SaveSolution(t, sol);
            SP.fullsave_t -= SP.fullsol_save_interval;
            step_since_dump = 0;
        }
        
        // Interface
        interface->CheckAll(t);
        if (interface->StopSimulationQ()){
            step_since_TV = -1;
            step_since_dump = -1; // Make sure it does final saves
        }
        
        
    }
    // Final energy and mean fields
    if (step_since_TV != 0){
        fluidEqs->Calc_Energy_AM_Diss(time_vars, t, sol);
        time_vars->Save_Mean_Fields(sol,  fft);
    }
    if (step_since_dump != 0){
        std::stringstream save_str;
        save_str << "Saving full solution at time " << t << std::endl;
        mpi.print1( save_str.str() );
        sol_save->SaveSolution(t, sol);
    }
    
    // Timing
    
    std::stringstream time_str;
    time_str << "Finished calculation at  t = " << t << " with average time-step " << integrator->mean_time_step() << std::endl << "Full time "<< interface->ElapsedTime() << "s:" << " IO time " << sol_save->IOtime() << "s: " << "Time variables: " << time_vars->TVtime() << "s" << std::endl;
    mpi.print1( time_str.str() );
    

    
    delete interface;
    delete fluidEqs;
    delete sol;
    delete integrator;
    delete time_vars;
    delete sol_save;
#ifdef USE_MPI_FLAG
    MPI_Finalize();
#endif
    return 0;
}

