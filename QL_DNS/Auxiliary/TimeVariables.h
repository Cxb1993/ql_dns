//
//  TimeVariables.h
//  QL_DNS
//
//  Created by Jonathan Squire on 8/25/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef __QL_DNS__TimeVariables__
#define __QL_DNS__TimeVariables__

#include "../General_Definitions.h"

#include "Input_parameters.h"
#include "fftwPlans.h"


//  Timevar storage, i.e., enregy, angular momentum etc.
class TimeVariables {
public:
    // Constructor
    TimeVariables(Inputs& SP,  int ENwidth, int numMF, int num_reynolds,int mpinode); // nsteps is number of saves, width is number to save
    // Destructor
    ~TimeVariables();
    
    // Increase time save number (could improve later)
    // Returns pointer to current array
    double* current_energy() {return energy_;};
    double* current_AM() {return angular_momentum_;};
    double* current_diss() {return dissipation_;};
    double* current_reynolds() {return reynolds_;};
    
    // Variables to calculate
    bool energy_save_Q() {return en_save_Q_;};
    bool AngMom_save_Q() {return AM_save_Q_;};
    bool dissip_save_Q() {return diss_save_Q_;};
    bool reynolds_save_Q() {return rey_save_Q_;};
    bool MeanField_save_Q() {return MF_save_Q_;};
    
    
    int number_of_saves() {return nsteps_;};
    
    // Save data - Change to each time-step?
    void Save_Data(double t);
    
    // Save mean fields - could improve
    void Save_Mean_Fields(dcmplxVec *MFdata,  fftwPlans& fft);
    
    
private:
    // Basic data
    const int nsteps_;  // Number of saves (length of array)
    const int width_;  // Number of saves for each step in each
    const int numMF_;   // Number of mean fields
    // Reynolds stress
    int num_reynolds_saves_; // Number of saves in Reynolds stress
    
    // Saving data booleans
    bool en_save_Q_,AM_save_Q_,diss_save_Q_;
    bool rey_save_Q_;
    bool MF_save_Q_;
    
    // Actual storage arrays
    double * energy_;
    double * angular_momentum_;
    double * dissipation_;
    double * reynolds_;
    // Saving these variables - file streams
    std::ofstream fileS_energy_;
    std::ofstream fileS_angular_momentum_;
    std::ofstream fileS_dissipation_;
    std::ofstream fileS_reynolds_;
    std::ofstream fileS_time_;
    
    // Saving these variables - file name
    std::string fname_energy_;
    std::string fname_angular_momentum_;
    std::string fname_dissipation_;
    std::string fname_reynolds_;
    std::string fname_time_;
    
    // Mean fields
    std::ofstream fileS_mean_fields_;
    std::string fname_mean_fields_;
    dcmplxVec MFdata_c_; // Space for taking fft and converting to real
    doubVec MFdata_d_;
    long MFlen_; // Set to zero once and change at first save
    
    // Directory
    std::string simulation_dir_;
    
    // Current position in the time dimension
    static int curr_pos_;
    
    
    // MPI node on which the data is stored
    int mpi_node_;
    bool mpi_node_bool_;// Boolean, 1 on root, 0 otherwise
    
};






#endif /* defined(__QL_DNS__TimeVariables__) */
