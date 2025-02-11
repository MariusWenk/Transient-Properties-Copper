/* 
 * File:   used_data.h
 * Author: Wenk, Ag Rethfeld
 *
 * Created on 16. May 2022, 11:21
 */

#ifndef DOS_H
#define DOS_H

#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <vector>
#include <string>
#include <functional>
#include <boost/math/tools/roots.hpp>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#include "integration.h"
#include "atom_unit_physics.h"
#include "physical_constants.h"
#include "used_data.h"

typedef std::vector<double> dv;

class DOS {

	public:
        DOS();
    	DOS(std::string dos_filename);
    	virtual ~DOS();

    	void read_DOS(std::string dos_filename);

        std::vector<std::vector<double>> dosfile; // Matrix or vector of vectors which will contain all 2 columns (E, sp+d_dos)

    //private:

    	//GSL accelerator of tot-dos
    	gsl_interp_accel *total_dos_accelerator;
    	//GSL spline of tot_dos
    	gsl_spline *total_dos_spline;
    	//GSL accel of sp_dos
   	 	gsl_interp_accel *sp_accel;
    	//GSL spline of sp_dos
    	gsl_spline *sp_spline;
    	//GSL accel fro d_dos
    	gsl_interp_accel *d_accel;
	    //GSL spline of d_dos
    	gsl_spline *d_spline;
		//Array for energy
    	double* energy;
    	//Array for total DOS
    	double* tot_DOS;
	    //Array for sp-DOS
    	double *sp_DOS;
    	//Array for d-DOS
    	double *d_DOS; 
    	double fermi_energy;
    	Cintegration *integrator_ptr;
    	double LowerBound;
   		double UpperBound;
		unsigned int n; // number of disc points
		double sp_band_mass;
		double d_band_mass;
 		double temperature_start;
		//! Distance from d-band edge to the Fermi-Energy
		double distance_d_band_E_F;
		double unitVolume;

};

#endif