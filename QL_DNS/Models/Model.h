//
//  Model.h
//  QL_DNS
//
//  Created by Jonathan Squire on 4/25/14.
//  Copyright (c) 2014 J Squire. All rights reserved.
//

#ifndef QL_DNS_Model_h
#define QL_DNS_Model_h

#include "../General_Definitions.h"
#include "../Auxiliary/MPIdata.h"
#include "../Auxiliary/TimeVariables.h"

//#include "../Auxiliary/TimeVariables.h"

// Abstract model class - container for equations of motion for QL_DNS
// Solution is stored in "solution" class

class solution;
class Kdata;

class Model {
public:
    
    Model(int NZ,const int NXY[2], const double * L) :
    L_(L)
    {
        // These dimensions do not depend on model (i.e., number of variables, MFs etc.)
        NZ_ = NZ;
        Nxy_[0] = NXY[0]; // Nyquist (and zero) frequncy included in NXY[0]
        Nxy_[1]= NXY[1]/2; // ONLY USE HALF OF THE TRANSFORM FOR THE Y DIMENSION DUE TO FULL Ckl BEING REAL!
        nz_dealiased_ = 1+2*int( NZ_/3.0 );
        
        nxy_full_= (Nxy_[0]-1)*Nxy_[1] ;// Leave out Nyquist frequency in kx
        // NB: Currently evolving the (0,0) component of Ckl. This is uneccessary but more convenient (especially with parallelization), will just stay zero anyway.
    }
    virtual ~Model() {};
    
    
    // fx = f(x,t)
    // Returns nonlinear part of of RHS and a linear operator (linop_Ckl), so
    //  as to allow for semi-implicit integrators.
    // The linear operators are evaluated at t + dt_lin
    virtual void rhs(const double t, const double dt_lin,
                    const solution * SolIn,solution * SolOut, doubVec **linOpFluct) = 0;
    
    // Initialization of linear operators - for mean fields linop is constant anyway
    virtual void linearOPs_Init(double t0, doubVec **linOpFluct, doubVec *linOpMF) = 0;
    
    // number of states (size of x)
    int Dimxy_full() const {return nxy_full_;}; // Full C size in x,y
    int NZ() const {return nz_dealiased_;};  // // Dealiased NZ dimension - this is dimension that is passed around with MPI
    int NZfull() const {return NZ_;}; // For speed/clarity
    virtual int Dimxy() const = 0; // Size of C array on given MPI process
    virtual int num_MFs() const = 0;  // Number of mean fields
    virtual int num_Lin() const =0; // Number of linear (fluctuating) fields
    
    // Box size
    double L(int ind) const { return L_[ind]; }; // Size of box
    int N(int ind) const { // Convenient to have full dimensions available
        if (ind==0){
            return Nxy_[0];
        } else if (ind==1){
            return Nxy_[1]*2; // Full y dimension in real space!
        } else if (ind==2){
            return NZ_;
        }
        return 0;
    };

    // Store k vectors as public member, only very rarely need to be accessed from outside model sub-class (for saving and restarting)
    Kdata *K;
//    
    // Remapping procedure
    virtual void ShearingBox_Remap(double qt, solution *sol) = 0;
//
//    // Quasi-linear state
//    virtual void set_QL_YN(bool QL) = 0;
//    
//
    ///////////////////////////////////////////////////////
    //  Generate Noise - since I'm not using a true stochastic integrator
    virtual void DrivingNoise(double t, double dt, solution *sol) = 0;
    ///////////////////////////////////////////////////////
    
    
    //////////////////////////////////////////////////////////////////
    //////  AUXILIARY FUNCTIONS OPERATING ON SOLUTION   //////////////
    virtual void Calc_Energy_AM_Diss(TimeVariables* tv, double t, const solution *sol) = 0;
    //////////////////////////////////////////////////////////////////
    
    ///////////////////////////////////////////////////
    // CFL condition - returns CFL number (e.g., kmax*By+q)
    virtual double Calculate_CFL(const solution *sol) = 0;
    ///////////////////////////////////////////////////

protected:
    
    // DIMENSIONS - these do not depend on model
    int Nxy_[2];           // Dimensions (Nx, Ny/2)
    int NZ_;
    int nz_dealiased_; // Dealiased NZ
    int nxy_full_;      // Ckl Nxy[0]*Nxy[1] total dimension
    const double * L_;       // Box size (Lx,Ly,Lz)
        
};




#endif  // MODEL_H_
