//
//  TimeVariables.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 8/25/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "TimeVariables.h"
#include "../solution.h" // Have to include here to stop circular definition problems




/////////////////////////////////////////////////////
//                                                 //
//                 TimeVariables                   //
//                                                 //
/////////////////////////////////////////////////////

// Current position
int TimeVariables::curr_pos_ = 0;

// Constructor
TimeVariables::TimeVariables(Inputs& SP, int ENwidth, int numMF, int num_reynolds, int mpinode) :
nsteps_(SP.nsteps/SP.timevar_save_nsteps + 1), width_(ENwidth), num_reynolds_saves_(num_reynolds),
mpi_node_(mpinode),
simulation_dir_(SP.simulation_dir), numMF_(numMF),nz_full_(SP.NZ),
en_save_Q_(SP.energy_save_Q), AM_save_Q_(SP.AM_save_Q), diss_save_Q_(SP.dissipation_save_Q),
rey_save_Q_(SP.reynolds_save_Q),
MF_save_Q_(SP.mean_field_save_Q),
clk_diff_(0.0)// Timing
{
    ////////////////////////////////
    //    THIS NEEDS TIDYING!!!!!!!!
    
    mpi_node_bool_ =  (mpi_node_ ==0);
    // Main variables
    // Store only on 1-processor - presumbably this is better
    // TODO: Only allocate necessary memory?
    if (mpi_node_bool_) {
        energy_ = new double[width_];
        angular_momentum_ = new double[width_];
        dissipation_ = new double[width_];
        reynolds_ = new double[num_reynolds_saves_];
    }
    // For mean fields
    MFdata_c_ = dcmplxVec::Zero( nz_full_ );
    MFdata_d_ = doubVec::Zero(nz_full_ );
    
    // Saving to disk - TODO improve to save at regular intervals?
    fname_energy_ = simulation_dir_ + "energy.dat";
    fname_angular_momentum_ = simulation_dir_ + "angular_momentum.dat";
    fname_dissipation_ = simulation_dir_ + "dissipation.dat";
    // Reynolds stress
    fname_reynolds_ = simulation_dir_ + "reynolds_stress.dat";
    // Mean fields
    fname_mean_fields_ = simulation_dir_ + "mean_fields.dat";
    // Time
    fname_time_ = simulation_dir_ + "time_vec.dat";
    
    
    std::ios_base::openmode o_mode;
    if (SP.start_from_saved_Q) {
        o_mode = std::ios::binary | std::ios::app;
    } else {
        o_mode = std::ios::binary | std::ios::trunc;
    }
    
    // REMOVE APPEND - OTHERWISE ITS ANNOYING!
    if (mpi_node_bool_) {
        
        
        if (en_save_Q_) {
            fileS_energy_.open(fname_energy_.c_str(), o_mode);
            if (!fileS_energy_.is_open() ) {
                std::cout << "WARNING: " << fname_energy_ <<  " file unable to be opened" << std::endl;
            }
        }
        if (AM_save_Q_) {
            fileS_angular_momentum_.open(fname_angular_momentum_.c_str(), o_mode);
            if (!fileS_angular_momentum_.is_open() ) {
                std::cout << "WARNING: " << fname_angular_momentum_ <<  " file unable to be opened" << std::endl;
            }
        }
        if (diss_save_Q_) {
            fileS_dissipation_.open(fname_dissipation_.c_str(), o_mode);
            if (!fileS_dissipation_.is_open() ) {
                std::cout << "WARNING: " << fname_dissipation_ <<  " file unable to be opened" << std::endl;
            }
        }
        // Reynolds stress
        if (rey_save_Q_) {
            fileS_reynolds_.open(fname_reynolds_.c_str(), o_mode);
            if (!fileS_reynolds_.is_open() ) {
                std::cout << "WARNING: " << fname_reynolds_ <<  " file unable to be opened" << std::endl;
            }
        }
        // Mean fields
        if (MF_save_Q_) {
            fileS_mean_fields_.open(fname_mean_fields_.c_str(), o_mode);
            if (!fileS_mean_fields_.is_open() ) {
                std::cout << "WARNING: " << fname_mean_fields_ <<  " file unable to be opened" << std::endl;
            }
        }
        
        //        Time data
        fileS_time_.open(fname_time_.c_str(), o_mode);
        if (!fileS_time_.is_open() ) {
            std::cout << "WARNING: " << fname_time_ <<  " file unable to be opened" << std::endl;
        }
        
        ////////////////////////////////////////
    }
    
    
}

// Destructor
TimeVariables::~TimeVariables() {
    if (mpi_node_bool_) {
        
        delete[] energy_;
        delete[] angular_momentum_;
        delete[] dissipation_;
        delete[] reynolds_;
    }
    
    
    fileS_energy_.close();
    fileS_angular_momentum_.close();
    fileS_dissipation_.close();
    fileS_reynolds_.close();
    fileS_mean_fields_.close();
    fileS_time_.close();

    
}


// Save the data to disk
void TimeVariables::Save_Data(double t){
    // Timing put outside of this function (in Calc_Energy_AM_Diss) using start_timing() and finish_timing()
    
    // Run only on root process
    if (mpi_node_bool_){
        
        if (en_save_Q_ ) {
            // Must write each individually since don't know that it is continuous in memory
            fileS_energy_.write( (char*) energy_, sizeof(energy_)*width_);
            // Read into matlab as a double (4*nsteps) array
        }
        if (AM_save_Q_ ) {
            // Must write each individually since don't know that it is continuous in memory
            fileS_angular_momentum_.write( (char*) angular_momentum_, sizeof(angular_momentum_)*width_);
            // Read into matlab as a double (4*nsteps) array
        }
        if (diss_save_Q_) {
            // Must write each individually since don't know that it is continuous in memory
            fileS_dissipation_.write( (char*) dissipation_, sizeof(dissipation_)*width_);
            // Read into matlab as a double (4*nsteps) array
        }
        if (rey_save_Q_ ) {
            // Must write each individually since don't know that it is continuous in memory
            fileS_reynolds_.write( (char*) reynolds_, sizeof(reynolds_)*num_reynolds_saves_);
            
            // Read into matlab as a double (num_reynolds_saves_*nsteps) array
        }
        
        // Write
        fileS_time_.write( (char*) &t, sizeof(double) );
        
    }

    
}


// Save Mean field data
void TimeVariables::Save_Mean_Fields(solution *sol, fftwPlans& fft){
    // Start timing
    clk_start_ = clock();
    
    // Could certainly be improved - just dumps the data into fname_mean_fields
    // IN MATLAB
    // Saves MF[i] sequentially from 0 to numMF. Data is double
    
    if ( MF_save_Q_ && mpi_node_bool_ ){
 
        // Save each MF variable sequentially
        for (int i=0; i<numMF_; ++i) {
            // Take fft - could be done faster with c_to_r fft
            fft.inverse(sol->pMF(i), &MFdata_c_);
            MFdata_d_ = fft.fac1D()*MFdata_c_.real();
            //Save
            fileS_mean_fields_.write( (char*) MFdata_d_.data(), sizeof(double)*nz_full_);
        }
    }
    
    // Timing
    clk_diff_ += (clock() - clk_start_ )/ (double)CLOCKS_PER_SEC;
}
//                                                 //
/////////////////////////////////////////////////////
