//
//  Input_parameters.h
//  QL_DNS
//
//  Created by Jonathan Squire on 8/23/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __QL_DNS__Input_parameters__
#define __QL_DNS__Input_parameters__

#include "../General_Definitions.h"

#include "../Auxiliary/MPIdata.h"

// Class for handling the input parameters and storing the data.
// This should be a convenient way to input data and have the methods required to process this in a basic way

// Really more like a struct

class Inputs {
public:
    // Constructor
    Inputs(const MPIdata& mpi,const std::string& input_file_name);
    
    
    
    ///////////////////////////////////////////////
    /////       USER INPUTS (through file)   //////
    
    // Grid size
    int NXY[2]; // X and Y grid
    int NZ ; // Z grid
    // Box dimensions
    double L[3];
    
    // Time domain
    double dt; // Might be ignored (for a variable time-step integrator)
    double CFL; // CFL number, might be ignored (for a fixed step integrator)
    double t_final; // Final time
    double t_initial; // Initial time
    double t_start; // Start time (different when loading from file)
    int i_start; // Start step (different from 0 when loading from file)
    
    // Saving
    double timvar_save_interval;   // Interval for saving energy, AM etc.
    double fullsol_save_interval;   // Interval for saving full solution
    
    // Dissipation
    double nu;   // Viscosity
    double eta;   // Resisitivity - could ignore for hydrodynamic
    // Shear
    double q;  // Shear rate
    double omega; // Rotation rate
    // Noise
    double f_noise;  // Driving noise (un-normalized)
    double noise_range_low; // Low driving noise cutoff (default 0)
    double noise_range_high; // High driving noise cutoff (default inf)
    bool drive_only_velocityQ; // Drive only velocity fluctuations
    bool drive_only_magneticQ; // Drive only magnetic
    // Shearing box
    bool remapQ;   // Whether to remap or not
    // Include quasi-linear feedback
    bool QuasiLinearQ;
    
    // Time variables to save
    bool energy_save_Q;
    bool AM_save_Q;
    bool dissipation_save_Q;
    bool reynolds_save_Q;
    bool mean_field_save_Q;
    
    // Initial conditions
    double initial_By;
    
    // Start using saved data
    bool start_from_saved_Q;
    
    // Check that equations are correct in main.cc
    std::string equations_to_use;
    
    
    ///////////////////////////////////////////////
    
    
    ///////////////////////////////////////////////
    /////      Derived quantities          ////////
    int nsteps; // Total number of steps
    int timevar_save_nsteps; // Steps before saving timevar
    int fullsol_save_nsteps; // Steps before full save
    bool fullsol_save_Q; // Whether or not there is a full save
    // Times to use in saving
    double fullsave_t;
    double timevar_t;

    
    // Directory for data storage
    std::string simulation_dir;
    
    
    
    
    
    
    ////////////////////////////////////////////////
    //////              FUNCTIONS             //////
    
    void initialize_();
    // Read from text file
    template <class T>
    T Read_From_Input_File_(const std::string& varstr, const std::string& full_file, T default_value);
    // Find input file .DSSinput
    std::string Find_Input_File(const std::string& directory,const std::string& input_name);
    
private:
    int mpi_node_; // Mpi ode for internal usage
    
};


#endif /* defined(__QL_DNS__Input_parameters__) */
