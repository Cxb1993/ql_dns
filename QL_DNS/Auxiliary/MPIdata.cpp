//
//  MPIfunctions.cc
//  MRIDSS
//
//  Created by Jonathan Squire on 5/7/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#include "MPIdata.h"

// Class for handling MPI related data


// Splits up the XY grid among processors
// Asigns to private variables myn_xyi_min_, myn_xyi_max_, which contain the nxy indices between which each processor should calculate
// e.g. Nx=16, Ny=16, 8*(16-1)=120 total linear fields. Split this among 8 processors as (0,14)(15,29)...
void MPIdata::Split_NXY_Grid(int nxy_full){
#ifdef USE_MPI_FLAG
    // Check that grid can be evenly divided among processors
    // If not, should use fewer processors!!
    if (my_node_==0) {
        if (nxy_full%total_nodes_ != 0) {
            std::cout << "<<<<< Error >>>>>" << std::endl <<
            "Number of MPI processes must be a multiple of the grid size!!" << std::endl;
            std::cout << "Grid size: " << nxy_full << ", Processors: " << total_nodes_ << std::endl;
            ABORT;
        }
    }
    
    nxy_per_node_ = nxy_full/total_nodes_; // Number per process
    
    myn_xyi_min_ = my_node_*nxy_per_node_;
    myn_xyi_max_ = (my_node_+1)*nxy_per_node_-1;
    
    std::cout << "Proc " << my_n_v() << ": " << minxy_i() << " to " << maxxy_i() << std::endl;
#else
    nxy_per_node_ = nxy_full;
    myn_xyi_min_ = 0;
    myn_xyi_max_ = nxy_per_node_ -1;
#endif
}



/////////////////////////////////////////////////
////        PRINTING
// NB: Have to pass by value if you want to use stringstream
void MPIdata::print1(std::string instr) const {
    // Print single statement from processor 0
    if (my_node_ == 0) {
        std::cout << instr;
    }
}

void MPIdata::printAll(std::string instr) const{
    // Prints processor number and instr out sequentially
    // Guaranteed to print each seperately so don't get overlapping output
    for(int i = 0; i < total_n_v(); i++) {
        Barrier();
        if (i == my_n_v()) {
            std::cout << "Proc " << my_n_v() << std::endl;
            std::cout << instr << std::endl;
        }
    }

    
}


//////////////////////////////////////////////////
//  MPI BARRIER
void MPIdata::Barrier() const {
#ifdef USE_MPI_FLAG
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}


//////////////////////////////////////////////////
/////    REDUCE FUNCTIONS (WRAPPERS FOR MPI_Reduce, MPI_AllReduce)
void MPIdata::SumReduce_doub(double* in_p, double* out_p, int size) {
    // MPI_Reduce wrapper for data of type double - convenient for multiple reasons
    // Always reduces to processor 0 right now
#ifdef USE_MPI_FLAG
    MPI_Reduce(in_p, out_p, size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    for (int i=0; i<size; ++i) {
        out_p[i] = in_p[i];
    }
#endif

}

void MPIdata::SumReduce_IP_doub(double* in_p, int size) {
    // MPI_Reduce wrapper for data of type double - convenient for multiple reasons
    // Always reduces to processor 0 right now
#ifdef USE_MPI_FLAG
    if (my_n_v() == 0) {
        MPI::COMM_WORLD.Reduce(MPI::IN_PLACE, in_p, size, MPI::DOUBLE, MPI::SUM, 0);
    } else {
        MPI::COMM_WORLD.Reduce(in_p, NULL, size, MPI::DOUBLE, MPI::SUM, 0);
    }
    
#endif
    
}

// Sum AllReduce double
void MPIdata::SumAllReduce_IP_double(double* in_p, int size) {
    // MPI_AllReduce wrapper for data of type double - convenient for multiple reasons
    // This one is in-place
#ifdef USE_MPI_FLAG
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, in_p, size, MPI::DOUBLE, MPI::SUM);
#endif
}

// Sum AllReduce double
void MPIdata::SumAllReduce_double(double* in_p, double *out_p, int size) {
    // MPI_AllReduce wrapper for data of type double - convenient for multiple reasons
    // This one is in-place
#ifdef USE_MPI_FLAG
    MPI::COMM_WORLD.Allreduce( in_p, out_p, size, MPI::DOUBLE, MPI::SUM);
#else
    for (int i=0; i<size; ++i) {
        out_p[i] = in_p[i];
    }
#endif
}

