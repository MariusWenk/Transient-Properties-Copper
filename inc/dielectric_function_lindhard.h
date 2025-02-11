/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   dielectric_function_lindhard.h
 * Author: ndione, Ag Rethfeld
 *
 * Created on February 27, 2017, 4:18 PM
 */

#ifndef DIELECTRIC_FUNCTION_LINDHARD_H
#define DIELECTRIC_FUNCTION_LINDHARD_H

#include <gsl/gsl_integration.h>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <stdlib.h> //Needed for exit()
#include <cstdlib>
#include <complex>

#include "dielectric_function.h"
#include "physical_constants.h"
#include "integration.h"
#include "gsl_function_pp.h"
#include "atom_unit_physics.h"
#include "fermi_distribution.h"




class dielectric_function_lindhard :public dielectric_function{
public:

    dielectric_function_lindhard();
    virtual ~dielectric_function_lindhard();

    // Redefine functions as they are virtual members of the base class dielectric_function


    double CalcRealPartEquilibrium(double omega);
    double CalcRealPartEquilibrium(double omega, double wavevector);
 
    

    double CalcImaginaryPartEquilibrium(double omega);
    double CalcImaginaryPartEquilibrium(double omega, double wavevector);

//    double GetIntegralOpticalLimit();
 void do_timestep();
   
private:
    
//    double LindhardImaginaryPart(double omega, double q_wavevector);

    //!Pointer to the electron band for properties distribution and band structure
    //Celectrons* band;
  
    //!Integration class pointer
    Cintegration* integrator;

    //!Integral that has to be calculated for the optical limit
    double integral_optical_limit;
    
    

};

#endif /* DIELECTRIC_FUNCTION_LINDHARD_H */

