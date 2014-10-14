//
//  InterfaceOutput.h
//  QL_DNS
//
//  Created by Jonathan Squire on 10/13/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __QL_DNS__InterfaceOutput__
#define __QL_DNS__InterfaceOutput__

#include "../General_Definitions.h"

// Store pointers to lots of useful bits
#include "MPIdata.h"
#include "FullSave_Load.h"
#include "TimeVariables.h"
#include "../Integrators/Integrator.h"
#include "Input_parameters.h"
#include "../solution.h"

//  Interface to output or save data on current status of the run when it finds a file exists, like in snoopy
// Searches in CURRENT_BASE_DIR
//  Also is a full simulation timer

// Currently only outputs time details

class InterfaceOutput {
public:
    InterfaceOutput(MPIdata *mpi, TimeVariables *tv, FullSave_Load *fsaves, Integrator *integ, Inputs *sp);
    ~InterfaceOutput(){};
    
    // Current time
    double ElapsedTime() {return (clock() - start_time_ ) / (double)CLOCKS_PER_SEC; };
    
    // Check for files
    void CheckAll(double t); // Output 1 to continue simulation
    
    // Stop simulation?
    bool StopSimulationQ(){ return stop_simulation_;};
    
private:
    // Stop simulation
    bool stop_simulation_;
    
    // Node for printing - since node zero does a bunch of stuff by itself, use another (except if only one is being used!)
    int mpi_node_for_print_;
    // Number of integrator steps between checking
    const int loop_check_interval_;
    int check_i_; // Current i for checking
    
    // Filenames to search for
    std::string status_name_;
    std::string dumpstop_name_;
    
    // Timing
    clock_t start_time_ ;// Timing
    
    //   Output functions for each utility
    void output_status(double t);
    
    // Test if file exists - inline for speed
    bool file_exists (const std::string& name) {
        if (FILE *file = fopen(name.c_str(), "r")) {
            fclose(file);
            return true;
        } else {
            return false;
        }
    }
    
    // Auxilliary classes
    MPIdata *mpi_;
    FullSave_Load *fsave_;
    TimeVariables *tvs_;
    Integrator *integ_;
    Inputs *sp_;
};

#endif /* defined(__QL_DNS__InterfaceOutput__) */
