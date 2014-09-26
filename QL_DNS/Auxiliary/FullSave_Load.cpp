//
//  FullSave_Load.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 9/15/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#include "FullSave_Load.h"

FullSave_Load::FullSave_Load(Inputs& SP, MPIdata* mpi, Model *eqs) :
saveQ_(SP.fullsol_save_Q),start_from_saved_Q_(SP.start_from_saved_Q),
mpi_(mpi),
K_(eqs->K),
clk_diff_(0.0)// Timing
{   // Sizes of arrays
    nz_ = eqs->NZ(); // This is dealiased size!!! (No point saving zeros)
    nxy_ = eqs->Dimxy();
    nxy_full_ = eqs->Dimxy_full();
    nlin_ = eqs->num_Lin();
    nMF_ = eqs->num_MFs();
    
    q_ = SP.q;
    
    // If specified that save is required, create file
    if (saveQ_ || start_from_saved_Q_) {
        // File name
        file_name_ = SP.simulation_dir + "FullSolution.h5";
        
        // Have to explicitly account for parallelism here
#ifdef USE_MPI_FLAG
        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
    
#ifdef USE_MPI_FLAG
        if (start_from_saved_Q_) {
            file = H5Fopen(file_name_.c_str(), H5F_ACC_RDWR, plist_id);
            if (file<0) {
                mpi_->print1("Unable to open .h5 file!!\n This is probably  caused be specifying start_from_saved_Q when no previous .h5 file exists.\n Proceding as per evaluation specified in input.");
                file = H5Fcreate(file_name_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
                start_from_saved_Q_ = 0;
            }
            
        } else {
            file = H5Fcreate(file_name_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
            if (file<0) {
                mpi_->print1("Unable to create .h5 file for saving!!");
            }
        }
        

#else
        if (start_from_saved_Q_) {
            file = H5Fopen(file_name_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
            if (file<0) {
                mpi_->print1("Unable to open .h5 file!!\n This is probably  caused be specifying start_from_saved_Q when no previous .h5 file exists.\n Proceding as per evaluation specified in input.");
                file = H5Fcreate(file_name_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
                start_from_saved_Q_ = 0;
            }

        } else {
            file = H5Fcreate(file_name_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            if (file<0) {
                mpi_->print1("Unable to create .h5 file for saving!!");
            }
        }
        
#endif
        // Buffers
        lin_write_buff_ = fcmplxVec(nlin_*nz_);
        MF_write_buff_ = fcmplxVec(nMF_*nz_);
        kbuff_ = new float[nxy_]; // Store real versions
        float *kybuff= new float[nxy_full_];
        if (mpi_->my_n_v() == 0){
            kxbuff_ = new float[nxy_full_];
        }
        
        // Sizes
        dimz_lin_ = new hsize_t[2];
        dimz_lin_[0] = nz_*nlin_;dimz_lin_[1] = 2;
        dimz_MF_ = new hsize_t[2];
        dimz_MF_[0] = nz_*nMF_;dimz_MF_[1] = 2;
        
        // Create datatype and space for the main saves
        // Fluctuating variables
        lin_dspace = H5Screate_simple(2, dimz_lin_, NULL); // Rank 1 with complex numbers -> 2
        lin_dtype = H5Tcopy(H5T_NATIVE_FLOAT); // Same for all arrays
        H5Tset_order(lin_dtype, H5T_ORDER_LE);
        lin_dset = new hid_t[nxy_full_];
        // Mean field variables
        MF_dspace = H5Screate_simple(2, dimz_MF_, NULL); // Rank 1 with complex numbers -> 2
        MF_dtype = H5Tcopy(H5T_NATIVE_FLOAT); // Same for all arrays
        H5Tset_order(MF_dtype, H5T_ORDER_LE);
        
        // kx stuff
        kx_dspace = H5Screate_simple(1, &nxy_full_, NULL);
        
        
        // Add ky and kz at this time. Could add other global properties (that don't change in time) here if this ends up being important
        htri_t kyexist = H5Lexists( file, "/ky", H5P_DEFAULT );
        htri_t kzexist = H5Lexists( file, "/kz", H5P_DEFAULT );
        // ky
        for (int j=0; j<nxy_; ++j)
            kbuff_[j] = K_->ky[j].imag();
        mpi_->PassToNode0_float(kbuff_, kybuff, nxy_); // Pass to 0 proc
        hid_t ky_dset;
        if (!kyexist) {
            ky_dset = H5Dcreate(file, "ky", H5T_NATIVE_FLOAT, kx_dspace,  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            SP.start_from_saved_Q = 0; // For TimeVariables - the only reason this wouldn't exist is if initial conditions are written by matlab, in which case should overwrite energy etc.
        } else {
            ky_dset = H5Dopen2(file, "/ky", H5P_DEFAULT);
        }
        
        // Write to file
        if (mpi_->my_n_v() == 0)
            H5Dwrite(ky_dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, kybuff );
        
        
        // kz
        hid_t kz_dspace  = H5Screate_simple(1, &nz_, NULL);
        hid_t kz_dset;
        if (!kzexist) {
            kz_dset = H5Dcreate(file, "kz", H5T_NATIVE_DOUBLE, kz_dspace,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        } else {
            kz_dset = H5Dopen2(file, "/kz", H5P_DEFAULT);
        }
        doubVec kzreal = (K_->kz).imag(); // Define real kz, much easier!
        H5Dwrite(kz_dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,kzreal.data() );
        
        
        // Close everything up
        H5Dclose(ky_dset);
        H5Dclose(kz_dset);
        H5Sclose(kz_dspace);
    
        delete[] kybuff;
    }
    
}


FullSave_Load::~FullSave_Load(){
    if (saveQ_ || start_from_saved_Q_) {
        
        H5Sclose(lin_dspace);
        H5Sclose(MF_dspace);
        H5Sclose(kx_dspace);
        H5Tclose(lin_dtype);
        H5Tclose(MF_dtype);
        
#ifdef USE_MPI_FLAG
        H5Fclose(file);
        H5Pclose(plist_id);
#else
        H5Fclose(file);
#endif
        
        delete[] dimz_lin_;
        delete[] dimz_MF_;
        delete[] kbuff_;
        delete[] lin_dset;
        if (mpi_->my_n_v() == 0){
            delete[] kxbuff_;
        }

    }
}


void FullSave_Load::SaveSolution(double t, solution *sol){
    // Dumps full solution into a new hdf5 group
    // Includes mean fields and fluctuations (indexed by kx,ky)
    
    // Start timing
    clk_start_ = clock();
    
    // Create group with time of save
    group_name_.str("");
    group_name_ << "/t=" << std::setprecision(15) << t;
    group_t = H5Gcreate2(file, group_name_.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // Create data sets on all processes
    ////////////////////////////////////////
    // Linear fields
    for (int jj=0; jj<nxy_full_; ++jj){ // jj symbolizes loop over global kx
        // data set name
        dsetname_.str("");
        dsetname_ << "Ind " <<jj <<": kx " << K_->kx_index_full[jj] << " ky " << K_->ky_index_full[jj];
        // Call HDF5 to create dataset
        lin_dset[jj] = H5Dcreate(group_t, dsetname_.str().c_str(), lin_dtype, lin_dspace,  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
    }
    
    //    Write Data
    // Fluctuating fields
    for (int i=0; i<nxy_; ++i) { // i symbolizes loop over local kx
        // Define full index
        int jj = i + mpi_->minxy_i();
        
        // Create buffer arrays
        for (int nV=0; nV<nlin_; ++nV) {
            lin_write_buff_.segment(nV*nz_, nz_) = (sol->pLin(i, nV))->cast<fcmplx>();
        }
        // Write to file
        H5Dwrite(lin_dset[jj], lin_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, lin_write_buff_.data() );
    }
    ////////////////////////////////////
    
    //////////////////////////////////////
    // Mean fields
    MF_dset = H5Dcreate(group_t, "Mean", MF_dtype, MF_dspace,  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // Create buffer array
    for (int nV=0; nV<nMF_; ++nV) {
        MF_write_buff_.segment(nV*nz_, nz_) = (sol->pMF(nV))->cast<fcmplx>();
    }
    // Write to file
    H5Dwrite(MF_dset, MF_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, MF_write_buff_.data() );
    //////////////////////////////////////

    //////////////////////////////////////
    // kx, ky array
    // This is done very badly! But, probably a small portion of total time (CHECK THIS!)
    // These copies could also be removed using derived mpi type!
    for (int j=0; j<nxy_; ++j)
        kbuff_[j] = K_->kx[j].imag()+q_*t*K_->ky[j].imag();
    mpi_->PassToNode0_float(kbuff_, kxbuff_, nxy_);
    // ky is stored in base of group (doesn't change with time)
    // Mean fields
    kx_dset = H5Dcreate(group_t, "kx", H5T_NATIVE_FLOAT, kx_dspace,  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Write to file
    if (mpi_->my_n_v() == 0)
        H5Dwrite(kx_dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, kxbuff_ );
    //////////////////////////////////////

    
    
    // Close datasets
    for (int jj=0; jj<nxy_full_; ++jj)
        H5Dclose(lin_dset[jj]);
    H5Dclose(MF_dset);
    H5Dclose(kx_dset);

    
    // CLose group
    H5Gclose(group_t);
    
    // Timing
    clk_diff_ += (clock() - clk_start_ )/ (double)CLOCKS_PER_SEC;
    
}



// Load solution from .h5 file
// Chooses last time-slice from .h5 file and loads it into the solution for a restart.
// Can also be used to load in initial conditions if necessary

// Resets Solution (initial condition), various parts of Inputs (starting time), and K->kx, to saved values
void FullSave_Load::LoadSolution_forRestart(solution *sol, Inputs &SP){
    
    // Start timing
    clk_start_ = clock();

    // Have already checked .h5 file exists in the constructor
    
    // Find the latest time (largest t) saved group
    std::string fsave_name;
    hid_t fsave_grp;
    double max_t = find_latest_time_save(file, &fsave_grp, fsave_name);
    fsave_name = "/"+fsave_name;
    mpi_->print1("Loading data for restart from " + fsave_name + "\n");
    SP.t_start = max_t;
    if (SP.fullsol_save_interval == SP.t_final-SP.t_initial) // Fix so that there's a final save
        SP.fullsol_save_interval = SP.t_final-SP.t_start - 1e-10;
    
    // Open data sets on all processes
    ////////////////////////////////////////
    // Linear fields
    for (int jj=0; jj<nxy_full_; ++jj){ // jj symbolizes loop over global kx
        // data set name
        dsetname_.str("");
        dsetname_ << fsave_name<< "/Ind " <<jj <<": kx " << K_->kx_index_full[jj] << " ky " << K_->ky_index_full[jj];
        // Call HDF5 to create dataset
        lin_dset[jj] = H5Dopen2(fsave_grp, dsetname_.str().c_str(), H5P_DEFAULT);
        
    }
    //    Read Data
    // Fluctuating fields
    for (int i=0; i<nxy_; ++i) { // i symbolizes loop over local kx
        // Define full index
        int jj = i + mpi_->minxy_i();
        
        // Read from file
        H5Dread(lin_dset[jj], lin_dtype, H5S_ALL, H5S_ALL,H5P_DEFAULT, lin_write_buff_.data());
        
        // Write into solution arrays
        for (int nV=0; nV<nlin_; ++nV) {
            *sol->pLin(i, nV)=lin_write_buff_.segment(nV*nz_, nz_).cast<dcmplx>();
        }
    }
    ////////////////////////////////////


    //////////////////////////////////////
    // Mean fields
    dsetname_.str("");
    dsetname_ << fsave_name<< "/Mean";
    MF_dset = H5Dopen2(fsave_grp, dsetname_.str().c_str(), H5P_DEFAULT);
    // Read from file
    H5Dread(MF_dset, MF_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, MF_write_buff_.data() );

    // Write into solution arrays
    for (int nV=0; nV<nMF_; ++nV) {
        *sol->pMF(nV) = MF_write_buff_.segment(nV*nz_, nz_).cast<dcmplx>() ;
    }
    mpi_->print1("Done loading solution\n");
    //////////////////////////////////////

    //////////////////////////////////////
    // kx array
    // This is done very badly! But, probably a small portion of total time (CHECK THIS!)
    dsetname_.str("");
    dsetname_ << fsave_name<< "/kx";
    htri_t kxexist = H5Lexists( fsave_grp, "kx", H5P_DEFAULT );
    if (kxexist) {
        kx_dset = H5Dopen2(fsave_grp, dsetname_.str().c_str(), H5P_DEFAULT);
        // Read from file
        if (mpi_->my_n_v() == 0)
            H5Dread(kx_dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, kxbuff_ );
        
        // Copy to other procs
        mpi_->ScatterFromNode0_float(kxbuff_, kbuff_, nxy_);
        // These copies could also be removed using derived mpi type!
        for (int j=0; j<nxy_; ++j)
             K_->kx[j] = dcmplx(0,kbuff_[j] - q_*SP.t_start*K_->ky[j].imag() );
        H5Dclose(kx_dset);
    } else {
        mpi_->print1("No kx data found, proceding using kx as initialized\n");
    }
    //////////////////////////////////////

    // Close up
    for (int jj=0; jj<nxy_full_; ++jj)
        H5Dclose(lin_dset[jj]);
    H5Dclose(MF_dset);
    
    
    H5Gclose(fsave_grp);
    
    // Timing
    clk_diff_ += (clock() - clk_start_ )/ (double)CLOCKS_PER_SEC;
    
}



// Find last time save in hdf5 file
// return t_max, time of last save
// file - file ID (should be opened in parallel if desired)
// fsave_grp - pointer to group id (output)
// fsave_name - name of group (output)
double FullSave_Load::find_latest_time_save(hid_t file, hid_t *fsave_grp, std::string &fsave_name){
    //  Base group
    hid_t basegroup = H5Gopen(file,"/",H5P_DEFAULT);
    // Number of objects in base
    hsize_t nobj;
	H5Gget_num_objs(basegroup, &nobj);
    
    // Return names of group
    char ** memb_names = new char*[(int)nobj];
    for (int i=0; i<nobj; ++i)
        memb_names[i] = new char[MAX_NAME];
    // Vectors to store index and id
    hid_t memb_ids[(int) nobj];
    std::vector<int> sub_grp_ind;
    
    // Loop through all objects in base
	for (int i = 0; i < nobj; i++) {
        // Find each object's name
		hsize_t len = H5Gget_objname_by_idx(basegroup, (hsize_t)i,
                                            memb_names[i], (size_t)MAX_NAME );
        
        // Find the type of the object
		int otype =  H5Gget_objtype_by_idx(basegroup, (size_t)i );
        if (otype == H5G_GROUP) {// Add to subgroup vector
            memb_ids[i] = H5Gopen(basegroup,memb_names[i],H5P_DEFAULT); // Non-group members will be uninitialized
            sub_grp_ind.push_back( i );
        }
    }
    
    // Loop through sub groups found in previous step - find that with the highest t value
    double max_t=-1.0;
    int max_i = -1 ;
    for (std::vector<int>::iterator it = sub_grp_ind.begin() ; it != sub_grp_ind.end(); ++it){
        // Name of group
        std::string gname(memb_names[*it]);
        // Position of t=, and number string
        double t;
        std::cin.precision( 16 );
        int pos = gname.find_first_of("t=");
        if (pos != 0) {
            mpi_->print1("Group name doesn't match t=... format\n");
        } else {
            std::string str_end = gname.substr(pos+2,pos+100);
            std::istringstream str_end_SS(str_end);
            str_end_SS >> t;
        }
        if (max_t < t) {
            max_t = t;
            max_i = *it;
        }
    }
    // Close groups all the other groups
    if (max_i >= 0) {
        for (std::vector<int>::iterator it = sub_grp_ind.begin() ; it != sub_grp_ind.end(); ++it){
            if (*it != max_i) {
                H5Gclose(memb_ids[*it]);
            }
        }
        *fsave_grp = memb_ids[max_i]; // Copy id to input pointer (don't copy pointer!)
        fsave_name = memb_names[max_i];
    } else {
        mpi_->print1("No time save found in .h5 file!");
        *fsave_grp = -1;
    }

    // Clean up
    for (int i=0; i<nobj; ++i)
        delete[] memb_names[i];
    delete[] memb_names;
    
    H5Gclose(basegroup);
    
    return max_t;
    
}

