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
    Barrier();
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
void MPIdata::SumAllReduce_doub(double* in_p, double* out_p, int size) {
    // MPI_Reduce wrapper for data of type double - convenient for multiple reasons
    // Always reduces to processor 0 right now
#ifdef USE_MPI_FLAG
    MPI_Allreduce(in_p, out_p, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    for (int i=0; i<size; ++i) {
        out_p[i] = in_p[i];
    }
#endif
    
}

void MPIdata::SumReduce_int(int* in_p, int* out_p, int size) {
    // MPI_Reduce wrapper for data of type int - convenient for multiple reasons
    // Always reduces to processor 0 right now
#ifdef USE_MPI_FLAG
    MPI_Reduce(in_p, out_p, size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    for (int i=0; i<size; ++i) {
        out_p[i] = in_p[i];
    }
#endif
    
}

void MPIdata::SumAllReduce_int(int* in_p, int* out_p, int size) {
    // MPI_Reduce wrapper for data of type int - convenient for multiple reasons
    // Always reduces to processor 0 right now
#ifdef USE_MPI_FLAG
    MPI_Allreduce(in_p, out_p, size, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
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


// Sum AllReduce dcmplx
void MPIdata::SumAllReduce_dcmplx(dcmplx* in_p, dcmplx *out_p, int size) {
    // MPI_AllReduce wrapper for data of type dcmplx - convenient for multiple reasons
#ifdef USE_MPI_FLAG
    MPI::COMM_WORLD.Allreduce( in_p, out_p, size, MPI::DOUBLE_COMPLEX, MPI::SUM);
#else
    for (int i=0; i<size; ++i) {
        out_p[i] = in_p[i];
    }
#endif
}

// Basic MPI send - complex data
void MPIdata::Send_dcmplx(dcmplx* s_buff, int count, int dest, int tag){
#ifdef USE_MPI_FLAG
    MPI_Send(s_buff, count, MPI::DOUBLE_COMPLEX, dest, tag, MPI_COMM_WORLD);
#else
    print1("SHOULD NOT BE CALLING Send_dcmplx WITHOUT MPI!!!\n");
#endif
}
// Basic MPI receive - complex data
void MPIdata::Recv_dcmplx(dcmplx* r_buff, int count, int source, int tag){
#ifdef USE_MPI_FLAG
    MPI_Recv(r_buff, count, MPI::DOUBLE_COMPLEX, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#else
    print1("SHOULD NOT BE CALLING Recv_dcmplx WITHOUT MPI!!!\n");
#endif
}



// Passes data from in_p to rec_p on processor 0
// rec_p must be of size np*size_init
void MPIdata::PassToNode0_float(float *in_p, float *rec_p, int size_init){
    // NB: This is essentially MPI_Gather, I hadn't realized and doesn't seem worth changing
#ifdef USE_MPI_FLAG
    for (int i=1; i<total_n_v(); ++i) {
        if (my_n_v() == i)
            MPI_Send(in_p, size_init, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
        if (my_n_v() == 0)
            MPI_Recv(rec_p+size_init*i, size_init, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
#endif
    // Copy over size_init elements from in_p to rec_p, same if no mpi
    if (my_n_v() == 0) {
        for (int i=0; i<size_init; ++i) {
            rec_p[i] = in_p[i];
        }
    }
}

// Scatters data from node 0 across all processes
void MPIdata::ScatterFromNode0_float(float *in_p, float *rec_p, int size_rec_p){
#ifdef USE_MPI_FLAG
    MPI_Scatter(in_p, size_rec_p, MPI_FLOAT,rec_p, size_rec_p*total_n_v(), MPI_FLOAT, 0,MPI_COMM_WORLD);
#else
    // Copy in_p to rec_p
    for (int i=0; i<size_rec_p; ++i) {
        rec_p[i] = in_p[i];
    }
#endif
}

void MPIdata::BroadcastFromNode0_doub(double *in_p, int size){
#ifdef USE_MPI_FLAG 
    MPI_Bcast(in_p, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}



