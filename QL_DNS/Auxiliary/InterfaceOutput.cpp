//
//  InterfaceOutput.cpp
//  QL_DNS
//
//  Created by Jonathan Squire on 10/13/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

//  Interface to output data about the current status of the run when it finds a file exists, like in snoopy
// Searches in CURR_BASE_DIR

// Note: DumpAndStop only works if full_save_ is not 0.0 in input file!!!! (Otherwise it hasn't done the necessary preparations for saving)

#include "InterfaceOutput.h"

InterfaceOutput::InterfaceOutput(MPIdata *mpi, TimeVariables *tv, FullSave_Load *fsaves, Integrator *integ, Inputs *sp) :
stop_simulation_(0),
mpi_node_for_print_(1),
loop_check_interval_(10),check_i_(0),// Step interval between checking for files
mpi_(mpi), fsave_(fsaves), tvs_(tv), integ_(integ), sp_(sp)
{
    // Timing
    start_time_ = clock();
    
    // File names to look for
    status_name_ = CURR_BASE_DIR + "status"; // Print status
    dumpstop_name_ = CURR_BASE_DIR + "DumpAndStop";  // Stop simulation and output save to h5 file (so it can be continued).
    
    // MPI - If only one node, change printing node to zero
    if (mpi_->total_n_v() == 1) {
        mpi_node_for_print_ = 0;
    }
    
}

//  Main function - check for files and do stuff if they're found
void InterfaceOutput::CheckAll(double t) {
    // Only chceck for file every loop_check_interval_ steps
    if (check_i_ == loop_check_interval_) {
        // Only check status on one node
        if ( mpi_->my_n_v() == mpi_node_for_print_ ) {
            //  Status
            if ( file_exists(status_name_) ){
                output_status(t);
            }
        }
        // Need to check for stop on all nodes
        if (file_exists(dumpstop_name_)) {
            // Stop simulation, saving current state (this is done automatically anyway)
            stop_simulation_ = 1;
            // Remove file
            mpi_->Barrier(); // Need barrier, since otherwise file might get deleted before other processes find it exists
            if (mpi_->my_n_v() == mpi_node_for_print_ ) {
                if( remove( dumpstop_name_.c_str() ) != 0 )
                    std::cout << "Error deleting status file!\n\n" ;
            }
        }
        
        check_i_ = 0 ;
    } else
        ++check_i_;

}


/////////////////////////////////////////////
//   These functions should only be executed on one process!
void InterfaceOutput::output_status(double t) {
    // Output info about current status of the run
    
    double diff = (clock() - start_time_ ) / (double)CLOCKS_PER_SEC;
    // Prints (to std::cout) info about the current run
    std::cout << "\nCalculation at t = " << t << " (final time t = " << sp_->t_final << ") after " << diff << "s." << std::endl << "Average time-step " << integ_->mean_time_step() << std::endl << "IO time " << fsave_->IOtime() << "s: " << "Time variables calculations: " << tvs_->TVtime() << "s\n" << "Current eta_K*k_max (U, B): (" << tvs_->etaK_times_kmax[0] << ", " << tvs_->etaK_times_kmax[1] << ")\n" << std::endl;
    
    // Delete file
    if( remove( status_name_.c_str() ) != 0 )
        std::cout << "Error deleting status file!\n\n" ;
    
}