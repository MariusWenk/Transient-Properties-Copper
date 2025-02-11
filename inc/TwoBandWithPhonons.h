/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TwoBandWithPhonons.h
 * Author: ndione
 *
 * Created on August 21, 2018, 3:17 PM
 */

#ifndef TWOBANDWITHPHONONS_H
#define TWOBANDWITHPHONONS_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h> 
#include <boost/numeric/odeint.hpp>

#include "atom_unit_physics.h"
#include "dielectric_function.h"
#include "dielectric_function_drude.h"
#include "fermi_distribution.h"
#include "gsl_function_pp.h"
#include "integration.h"
#include "interpolation1d.hpp"
#include "interpolation2d.hpp"
#include "Laser.h"
#include "physical_constants.h"
#include "root_finding.h"
#include "ode.h"
#include "solver.h"

typedef std::vector<double> state_type;
typedef unsigned int uint;


class TwoBandWithPhonons : public ODE{
    
	public:
		TwoBandWithPhonons();    
		
		TwoBandWithPhonons(std::string filename, Cfermi_distribution *fermi_ptr); //Use this constructor if the e-ph coupling input file is 1D

		//TwoBandWithPhonons(std::string sp_filename2D, std::string d_filename2D, Cfermi_distribution *fermi_ptr, uint col0, uint col1, uint col2); //Use this constructor if the e-ph coupling input file is 2D
   
	    	virtual ~TwoBandWithPhonons();

	    	std::vector<std::vector<double>> LoadData;
	    	//! Hand over the derivative of a physical quantity with respect to time
	    	void calculate_dxdt(double t);
	    	//! Hands over the starting time for ODEs            
	    	double get_t_start();
	    	//! Hands over the ending time for ODEs
	    	double get_t_end();
	    	//! Hands over the timesped  for ODEs
	    	double get_t_step();
	    	//!Hands over the number of photons used to excite the material
	    	double OpticalExcitation(double time);
	    	//!Hands over relaxtion terms of electrons
	    	double RelaxationTerm();
	    	//!Hands over the total absorbed energy from the laser
	    	double AbsorbedEnergy(double time);
	    	//! Hands over dynamnics  of conduction band density  
	    	void dn_sp_dt(double time);
	    	//! Hands over dynamics of d electrons , 
	    	void dn_d_dt(double time);
	    	//! Hands over the density of valence electrons (sp and d electrons), always return 0 
	    	void dn_val_dt(double time);
	    	//! Hands over internal energy of valence electrons
	    	void du_val_dt(double time);
	    	//! Hands over the varation of the chemical potential of sp electrons
	    	void dmu_sp_dt(double time);
	    	//! Hands over the varation of the chemical potential of  electrons
	    	void dmu_d_dt(double time);
	    	//! Hands over the varation of  equilibrium chemical potential  
	    	void dmu_eq_dt(double time);
	    	//!Hands over the dynamics of electrons temperature
	    	void dTe_dt(double time);
	    	//! Hands over change of the phonons temperature
	    	void dTph_dt(double time);
	    	//! Hands over the fermi energy of the material
	    	double GetFermiEnergy();
	    	//! Hands over the initial density of the conduction band
	    	double Get_n_sp_initial();
	    	//! Hands over the initial density of the conduction band
	    	double Get_n_d_initial();
	    	// ! Hands over the total density of the three bands sp, d and f
	    	double GetTotalDensity();
	    	//! Hands over relaxation time
	    	double GetRelaxationTime();
	    	//! Hands over the relaxation rate
	    	double GetRelaxRate();
	    	//! Hands over the unit cell volume of the material
	    	double GetUnitCellVolume();
	    	//! Hands over the total absoprtion coefficient
	    	double GetAlphaTot();
		//double CalcAlphaTot(double Te, double Tph);
		//double CalcAlphaTot_neq(double Tph, double sp_density, double d_density);
	    	//Hands over the percentage of interband absorption coefficient calculated a each time step
	    	double InterbandPercentageFromDynamics();
	    	//Hands over the percentage of intraband absorption coefficient calculated a each time step
	    	double IntrabandPercentageFromDynamics();
	    	//Hands over the starting temperature 
	    	double GetStartTemperature();
	    	//! Hands over starting temperature for the lattice
	    	double GetStartTemperature_ph();
	    	//! Hands over a temperature dependent e-ph-coupling by interpolating Lin's data
	    	double TeDependentElectronPhononCoupling(double electron_temperature);
	    	//! Hands over heat capacity of phonons
	    	double PhononHeatCapacity();
	    	//! Energy loss term. Electrons couple to the phonons and exchange energy.
	    	double GetLossTerm();
	 	//Return private variable for TPA
		double GetTwoPhotonAbsorptionCoeff();
		// Hands over the TPA coeff times laser intensity
		double TwoPhotonAbsorption(double time);
		double TwoPhotonAbsorptionFunction(double time);	     
		/*!Hands over a Te and density resolved coupling
		@param  number of electron
		@param  Te
		*/
		
		double impactIonization();
		double impactIonizationRate();
		double electronHoleRecombination();
		double electronHoleRecombinationRate();
    
	private:

		//Array to store electron temperature for interpolation
		double *Te_array;
		//Array to store electron-phonon coupling for interpolation
		double *e_ph_array;
		//! GSL acccelerator of e_ph
		gsl_interp_accel *e_ph_accel;
		//! GSL spline of e_ph
		gsl_spline *e_ph_spline;
		//!Pointer for properties of dielectric function and optical parameters
		dielectric_function *dielectric_ptr; // class pointer 
		Claser *laser_ptr; // Class pointer
		//!Pointer to the fermi distribution class
		Cfermi_distribution *fermi_ptr;
		//! Integration pointer
		Cintegration *integrator_ptr;
	
		 


	    	//! Fermi energy of material
	    	double fermi_energy;
	    	//!relaxation time
	    	double relax_time;
	    	//! Conduction electrons relaxation rate
	    	double relax_rate;
	    	//!unit_cell_volume
	    	double unit_cell_volume;
	    	//! total density of all 3 bands
	    	double total_density_2_bands;
	    	//! Initial sp density
	    	double n_sp_initial;
	    	//! Initial d density
	    	double n_d_initial;
	    	//! Total absorption coefficient from https://refractiveindex.info (Should be calculatede once equilibrium optical properties known)
	    	double alpha_tot;
	    	//! Interband absorption coefficient from https://refractiveindex.info
	    	//! starting temperature
	    	double elec_start_temperature;// for electrons
	    	double ph_start_temperature;
	
	    	double t_start;
	    	double t_end;
	    	double t_step;

		//Number of disc points of The DOS file
		int n;
		double twoPhotonAbsorptionCoeff; // TPA coeff .uses a constant for now
	 	ufd::Interpolation1D *interp_pointer1D;
		ufd::Interpolation2D *sp_interp_pointer2D;
		ufd::Interpolation2D *d_interp_pointer2D;

		uint col0, col1, col2;

		//Alias
		double& n_sp = x[0], & n_d = x[1], & n_val = x[2];
		double& u_val = x[3], & T_e = x[4], & T_ph = x[5];
		double& mu_sp = x[6], & mu_d = x[7], & mu_eq = x[8];
};

#endif /* TWOBANDWITHPHONONSH */


