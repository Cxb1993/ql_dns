//
//  MPIfunctions.h
//  MRIDSS
//
//  Created by Jonathan Squire on 5/7/14.
//  Copyright (c) 2014 Jonathan Squire. All rights reserved.
//

#ifndef __MRIDSS__MPIdata__
#define __MRIDSS__MPIdata__

#include "../General_Definitions.h"



class MPIdata {
public:
    // Need to create MPIdata object before knowing number of processors
    MPIdata(void) : my_node_(0), total_nodes_(1), communicator_size_(1) {}; // The my_node_ and total_nodes_ objects are changed by the call to MPI_Comm_Size etc. in main.cc!!!
    ~MPIdata(){};
    
    // Getter functions
    int* total_n_p() { return &total_nodes_; }; // Pointer to total_nodes
    int total_n_v() const { return total_nodes_; }; // Value of total number of nodes
    int* my_n_p() { return &my_node_; };  // Pointer to my_node
    int my_n_v() const { return my_node_; };    // Value of my node
    int* comm_size_p() { return &communicator_size_; }; //Pointer to communicator size whatever this is

    int nxy() const {return nxy_per_node_; };
    int minxy_i() const { return myn_xyi_min_; };
    int maxxy_i() const { return myn_xyi_max_; };
    
    
    // Set up splitting between nodes
    void Split_NXY_Grid(int nxy_full);
    
    // Printing
    void print1(std::string instr) const;
    void printAll(std::string instr) const;
    
    // Barrier
    void Barrier() const;
    
    
    ////////////////////////////////////////////////////
    // MPI reduce functions
    void SumReduce_doub(double* in_p, double* out_p, int size);// Sum double values
    void SumAllReduce_doub(double* in_p, double* out_p, int size);// Sum double values
    void SumReduce_int(int* in_p, int* out_p, int size); // Sum int values
    void SumAllReduce_int(int* in_p, int* out_p, int size); //Sum int values
    void SumReduce_IP_doub(double* in_p, int size);// Sum double values
    void SumAllReduce_dcmplx(dcmplx* in_p, dcmplx*out_p, int size);//Sum double,
    
    // Send and receive for noise
    void Send_dcmplx(dcmplx* s_buff, int count, int dest, int tag);
    void Recv_dcmplx(dcmplx* r_buff, int count, int source, int tag);
    // Pass all of an array to node 0
    void PassToNode0_float(float *in_p, float *rec_p, int size_init);
    // Scatter data from node 0
    void ScatterFromNode0_float(float *in_p, float *rec_p, int size_rec_p);
    // Broadcase from node zero
    void BroadcastFromNode0_doub(double *in_p, int size);
private:
    
    // General MPI data
    int total_nodes_;  // Number of processors
    int my_node_;     // My node
    int communicator_size_; // Something - might want later
    
    // Data for each node in nxy_full array
    int nxy_per_node_; // Number of CKl per node
    // Min/max index in i for nxy_full for each processor
    int myn_xyi_min_, myn_xyi_max_;
    
    
};

#endif /* defined(__MRIDSS__MPIdata__) */