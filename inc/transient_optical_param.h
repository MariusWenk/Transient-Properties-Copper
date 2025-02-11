 

/* 
 * File:     transient_optical_param.h
 * Author: Ndione. Ag Rethfeld
 * 
 * Created on 12. Juni 2018, 17:33
 */

#ifndef TRANSIENT_OPTICAL_PARAM_H
#define TRANSIENT_OPTICAL_PARAM_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h> 

#include "atom_unit_physics.h"
#include "dielectric_function.h"
#include "dielectric_function_drude.h"
#include "Laser.h"
#include "physical_constants.h"
#include "used_data.h"

class transient_optical_param {

	public:
		transient_optical_param();
    		transient_optical_param(std::string input_optics_file,   
    		dielectric_function_drude *drude_point, dielectric_function_lorentz *lorentz_point);
    		virtual ~transient_optical_param();

		//Calculate optics and write data to file
    		void OpticalPropertiesDrudeLorentz();
	    	//Hands over the laser frequency in eV
	    	double GetOmegaProbe();
	    	//!Hands over the laser wavelength
	    	double GetWavelength();
	    	//store time, n_sp, n_d...needed for optics
	    	std::vector<std::vector<double>> opt_param;

	private:
		//pointer to class drude_dielctric_function
		dielectric_function_drude *drude_point;
		//pointer to class lorentz_dielectric_function
		dielectric_function_lorentz *lorentz_point;
	    
	    	//probe frequency in eV
	    	double omega_probe;
	    	//probe wavelength in nm
	    	double wavelength;

	    	//Arrays to store time and quantities needed for optics
	    	double *Time;
	    	double *n_sp;
	    	double *n_d;
		double *T_lattice;
		double *T_electron;
		double *mu_sp;
		double *mu_d;
		double *mu_eq;
		
		//optics are saved in these arrays
		//complex dielctric functino
		std::complex<double>* complexDielectricFunction;
		//complex refractive index
		std::complex<double>* n_complex;
		//equilibrium dielectric function
		std::complex<double>* complexDielectricFunction_eq;
		double *omega;
		double *T_film_eq;
		double *R_film_eq;
		//real part of AC conductivity
		double *sigma_r;
		//imaginary part of AC conductivity
		double* sigma_i;
		//DC conductivity
		double* dcConductivity;
		//e-ph collision frequency
		double* e_ph_nu;
		//el-el collision frequency
	    	double* spd_nu;
	    	//total collision frequency
	    	double* tot_nu;
    	
	    	//Reflection, transmission, Absorption of bulk
	    	double* R_n_bulk;
	    	double* R_s_polarized_bulk;
	    	double* R_p_polarized_bulk;

		double* T_n_bulk;
		double* T_s_polarized_bulk;
		double* T_p_polarized_bulk;
		
		double* A_n_bulk;
		double* A_s_polarized_bulk;
		double* A_p_polarized_bulk;
		
		// Reflection, transmission, Absorption of a free standing thib film
		double* R_n_film;
	    	double* R_s_polarized_film;
	    	double* R_p_polarized_film;
	    	
	    	double* T_n_film;
		double* T_s_polarized_film;
		double* T_p_polarized_film;
		
		double* A_n_film;
		double* A_p_polarized_film;
		double* A_s_polarized_film;
		
		//size of all arrays
		unsigned int n_size;

		unsigned int probe_wavelength_amount;
		double min_probe_wavelength;
		double max_probe_wavelength;

		unsigned int time_step_factor;



};


#endif /* TRANSIENT_OPTICAL_PARAM_H */
