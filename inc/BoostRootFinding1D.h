/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BoostRootFinding1D.h
 * Author: Ndione, Ag Rethfeld
 *
 * Created on October 31, 2018, 1:52 PM
 */

#ifndef BOOSTROOTFINDING1D_H
#define BOOSTROOTFINDING1D_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <functional>
#include <boost/math/tools/roots.hpp>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#include "integration.h"
#include "physical_constants.h"
#include "atom_unit_physics.h"
#include "used_data.h"
 

class BoostRootFinding1D {

	public:

    		BoostRootFinding1D(std::string dos_filename);
		BoostRootFinding1D();
   		virtual ~BoostRootFinding1D();
    
    		std::vector<std::vector<double>> dosfile; // Matrix or vector of vectors which will contain all 2 columns (E, sp+d_dos)
    		double GetFermiDistribution(double energy, double T, double mu);
    		double NumberofParticles(double T, double mu);
    		double sp_NumberofParticles(double T, double mu);
    		double d_NumberofParticles(double T, double mu);
    		double df_dmu(double energy, double T, double mu);
    		double Get_P_mu(double T, double mu);
    		double GetMuRootFinding(double T);
    		double GetDOS(double energy);
		double GetDOS_sp(double energy);
    		double GetDOS_d(double energy);
    		double GetLowerBound();
    		double GetUpperBound();
		double GetFermiEnergy();
		double GetNumberDiscPoints();
		double sp_FEG(double energy);
		double d_FEG(double energy);
		double GetTemperatureStart();
		double GetDistance_d_BandtoFermi();
		double SmoothFermiFunction(double energy);
		double StepFunction(double energy);
		double GetUnitVolume();
		double get_sp_bandMass();
		double get_d_bandMass();

	private:

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
    		double* DOS;
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

#endif /* BOOSTROOTFINDING1D_H */


