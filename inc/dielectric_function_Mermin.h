/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   dielectric_function_Mermin.h
 * Author: ndione, Ag Rethfeld
 *
 * Created on 27. März 2018, 15:14
 */

#ifndef DIELECTRIC_FUNCTION_MERMIN_H
#define DIELECTRIC_FUNCTION_MERMIN_H


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


class dielectric_function_Mermin {
public:
    dielectric_function_Mermin();
    dielectric_function_Mermin(const dielectric_function_Mermin& orig);
    virtual ~dielectric_function_Mermin();
private:

};

#endif /* DIELECTRIC_FUNCTION_MERMIN_H */

