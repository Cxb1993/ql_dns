//
//  FullSave_Load.h
//  QL_DNS
//
//  Created by Jonathan Squire on 9/15/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __QL_DNS__FullSave_Load__
#define __QL_DNS__FullSave_Load__

#include "../General_Definitions.h"

#include "Input_parameters.h"
#include "MPIdata.h"
#include "../solution.h"
#include "../Models/Model.h"
#include "Kdata.h"

// Class to handle saving the full solution in hdf5 format. Also handles a complete restart of the simulation from some save


class FullSave_Load {
public:
    FullSave_Load(Inputs& SP, MPIdata* mpi, Model *eqs);
    ~FullSave_Load();
    
    // Store a time slice of the full data
    void SaveSolution(double t, solution *sol);
    
    // Load final time slice for a restart
    void LoadSolution_forRestart(solution *sol, Inputs &SP);
    
    bool start_from_saved(){return start_from_saved_Q_;}; // Whether to start from saved solution
    
    double IOtime() {return clk_diff_;}; // Accumulated time for full saves
private:
    const bool saveQ_;
    bool start_from_saved_Q_;
    
    std::string file_name_;
    std::stringstream group_name_;
    std::stringstream dsetname_;

    // Timing
    clock_t clk_start_;
    double clk_diff_; // Add up timing from all calls to check total is not dominating
    
    // Various hdf5 storage handles
    hid_t    plist_id;  // Parallel property
    hid_t       file;         /* file handle */
    hid_t       group_t;   /* group for t=0 etc. data set */
    hid_t       lin_dtype, lin_dspace; // data type and space for fluctuating variables
    hid_t   MF_dtype, MF_dspace;
    hid_t   *lin_dset; // data set (need to define all on all processors
    hid_t  MF_dset;// MF data
    hid_t kx_dspace, kx_dset;
    
    // Size storage
    hsize_t nz_, nxy_full_, nxy_;
    hsize_t nlin_, nMF_; // Number of mean and fluctuating fields
    hsize_t *dimz_lin_, *dimz_MF_; // Includes 2 for complex data
    
    // Pointers to other data structures
    MPIdata *mpi_;
    Kdata *K_;
    double q_;
    
    // Buffer for writing - not the fastest probably, but otherwise presumably need HDF5 appending capability which seems really annoying. Have to change type anyway, so may not make that much difference
    fcmplxVec lin_write_buff_;
    fcmplxVec MF_write_buff_;
    float *kxbuff_, *kbuff_; // Pass data to p0 to save, saves using hyperslabs but should upgrade!
    
    // Find last time save in hdf5 file
    double find_latest_time_save(hid_t file, hid_t *fsave_grp, std::string &fsave_name);
    
};
#endif /* defined(__QL_DNS__FullSave_Load__) */
