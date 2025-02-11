 

/**************************************** 
 * File:   FourBandWithPhonons.h	*
 * Author: Ndione, AG Rethfeld		*
 *					*
 * Created on May 06, 2020, 01:08 PM
 ***************************************/

#ifndef FOURBANDWITPHONONS_H
#define FOURBANDWITPHONONS_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
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

#include "atom_unit_physics.h"
#include "dielectric_function.h"
#include "fermi_distribution.h"
#include "gsl_function_pp.h"
#include "integration.h"
#include "interpolation1d.hpp"
#include "Laser.h"
#include "physical_constants.h"
#include "root_finding.h"
#include "ode.h"
#include "solver.h"
#include "spline.h"
#include "ThreeBandWithPhonons.h"


class FourBandWithPhonons :  public ODE  {

	public:
 
		FourBandWithPhonons();
 		FourBandWithPhonons(std::string filename, Cfermi_distribution *fermi_ptr);
     
    		virtual ~FourBandWithPhonons();
		
		std::vector<std::vector<double>> LoadData;
	
	    	//! Hand over the derivative of a physical quantity with respect to time
	    	void calculate_dxdt(double time);
	    	//! Hands over the starting time for ODEs            
	    	double get_t_start();
	    	//! Hands over the ending time for ODEs
	    	double get_t_end();
	    	//! Hands over the timesped  for ODEs
	    	double get_t_step();
	    
	    	//!Hands over the number of photons used to excite the material
	    	double OpticalExcitation(double time);
	    	//! Hands over relaxtion terms of electrons
	    	double RelaxationTerm();
	    	//! Hands over term accounting for the auger processs
	    	double AugerTerm_4fBand();
		double AugerTerm_5pBand();
	    	//! Hands over potential energy of valence electrons
	    	double PotentialEnergy(double time);
	    	//! Hands over Energy lost by valence electrons due to Auger recombination
	    	double AugerEnergy5p_4f();
	    	double AugerEnergy5d_5p();
	    	//!Hands over the total absorbed energy from the laser
	    	double AbsorbedEnergy(double time);
	    	//! Hands over dynamnics  of conduction band density  
	    	void dn_6sp_dt(double time);
	    	//! Hands over dynamics of d electrons 
	    	void dn_5d_dt(double time);
		//! Hands over dynamics of 5p electrons 
		void dn_5p_dt(double time);
	    	//! Hands over dynamics of f electrons 
	    	void dn_4f_dt(double time);
	    	//! Hands over the density of valence electrons (sp and d electrons), always return 0 
	    	void dn_val_dt(double time);
	    	//! Hands over internal energy of valence electrons
	    	void du_val_dt(double time);
	    	//! Hands over the varation of the chemical potential of sp electrons
	    	void dmu_6sp_dt(double time);
	    	//! Hands over the varation of the chemical potential of  electrons
	    	void dmu_5d_dt(double time);
	    	//! Hands over the varation of  equilibrium chemical potential  
	    	void dmu_eq_dt(double time);
	    	//!Hands over the dynamics of electrons temperature
	    	void dTe_dt(double time);
	    	//! Hands over the phonon temperature
	    	void dTph_dt(double time);
	    	//! Hands over the fermi energy of the material
	    	double GetFermiEnergy();
	    	//! Hands over the chemical potential
	    	double Get_n_6sp_initial();
	    	//! Hands over the initial density of the conduction band
	    	double Get_n_5d_initial();
	    	//! Hands over the initial density of the f shell
	    	double Get_n_4f_initial();
	  
	    	//! Hands over the the minimum energy of the conducton band
	    	double Get_E_6sp_min();
	    	//! Hands over the minimum energy of the valence band
	    	double Get_E_5d_min();
	    	//! Hands over the minimum energy of f shell
	    	double Get_E_4f_min();
		//! Hands over the minimum energy of f shell
	    	double Get_E_5p_min();
	    	//! Hands over the distance from the 4f shell up to 6sp
	    	double GetDistance6sp_4f();
		//! Hands over the distance from the 4f shell up to 5p
	    	double GetDistance6sp_5p();
	    	//! Hands over inherent linewidth of 4f electrons measurements
	    	double InherentLinewidth();
	    	//! Hands over lifetime calculated for the inherent linewidth
		double GetLifetime4fStates();
	    	//! Hands over the auger rate 
		double GetRate4fStates();
	    	//! Hands over relax lifetime
	    	double GetRelaxationLifetime();
	    	//! Hands over the relaxation rate
	    	double GetRelaxRate();
	    	//! Hands over the unit cell volume of the material
	    	double GetUnitCellVolume();
	    	//! Hands over the total absoprtion coefficient
	    	double GetAlphaTot();
	    	//! Hands over the interband absoprtion coefficient
	    	double GetAlphaInter();
	    	//! Hands over temperature at t=0
	    	double GetStartTemperaure();
	    	//! Hands over the loss term! this term is not fully realistic but may incorporate transport, phonons effeccts, radiation ...
    		double GetLossTerm();  
	    	//! Hands over starting temperature for the lattice
	    	double GetStartingPhononTemperature();
	    	//! Hands over a temperature dependent e-ph-coupling by interpolating Lin's data
	    	double TeDependentElectronPhononCoupling(double electron_temperature);
	    	//! Hands over heat capacity of phonons
	   	double PhononHeatCapacity();
		//!Hands over the electron density in th 5p-band
		double Get_n_5p_iniial();
		//!Hands lifetime of 5p-holes
		double GetLifetime5pStates();
		double GetRate5pStates();


	private:

		//Array to store electron temperature for interpolation
		double *Te_array;
		//Array to store electron-phonon coupling for interpolation
		double *e_ph_array;
		//! GSL acccelerator of e_ph
		gsl_interp_accel *e_ph_accel;
		//! GSL spline of e_ph
		gsl_spline *e_ph_spline;
	

	 	//!Pointer to the fermi distribution class
	    	Cfermi_distribution *fermi_ptr;
	    	//! Integration pointer
	    	Cintegration *integrator_ptr;
	    	//! Pointer for the laser
	    	Claser *laser_ptr;
	    	// Pointer to the dielectric class to access opt paramrs
	    	dielectric_function *dielectric_ptr;
	    	//! Fermi energy of material
	    	double fermi_energy;
	    	double start_temperature;  // for electrons 
	    	double ph_start_temperature;

	    	/*The lifetime broadening of 4f electrons in gold may be extracted using Heisenberg's uncertainty principle, as we know the experimental value of the 
	     	*inherent linewidths. The linewidth values include both the lifetime broadening and the fine structure contribution.
	     	* Experimental values of 4f binding energy and inherent linewidths in gold can be fing in this paper:
	     	*https://www.sciencedirect.com/science/article/pii/S0368204810000113
		*The linewidht of a level may be used to extract the lifetime of a single hole in that level using the Heisenberg uncertainty principle, namely \Gamma\tau=\hbar 
		*/
	    	double inherent_linewidth;
	    	//!Lifetime extracted from inherent linewidth 
	    	double lifetime_4f_states;;
	    	//! Three body recombination rate
		double rate_4f_states;
	    	//! Relaxation time
	    	double relax_time;
	    	//! Conduction electrons relaxation rate
	    	double relax_rate;
	    	//!unit_cell_volume
	    	double unit_cell_volume;

	    	//! Initial sp density
	    	double n_6sp_initial;
	    	//! Initial d density
	    	double n_5d_initial;
	    	//! Initial f density
	    	double n_4f_initial;         
		//!Initial density of 5p-band
		double n_5p_initial;
	    	//! Total absorption coefficient from https://refractiveindex.info 
	    	double alpha_tot;
	    	//! Interband absorption coefficient from https://refractiveindex.info
	    	double alpha_inter;
	    	//! factor obtain while rewritting density into chemical potential
	    	double mu_factor;

	    	//! Minimum energy conduction band
	    	double E_6sp_min;
	    	//! Minimum energy valence band
	    	double E_5d_min;
		//! Minimum energy 5p band
	    	double E_5p_min;	
		//! Distance between 6sp band and 5p shell
	    	double dist_6sp_5p;
	    	//! Distance between 6sp band and f shell
	    	double dist_6sp_4f;
	    	/* Minimum energy f shell
	     	*Core levels in XPS use the nomenclature nl_{j} where n is the principal quantum number, l is the angular momentum quantum number and j = l + s
	     	*All orbital levels except the s levels (l = 0) give rise to a doublet with the two possible states having different binding energies.
	     	* This is known as spin-orbit splitting. The peaks will also have specific area ratios based on the degeneracy of each spin state, 
	     	* i.e. the number of different spin combinations that can give rise to the total j. For our model, the 4f electrons are considered
	     	which leads to j = 3+1/2=7/2 and j=3-1/2=5/2. Hence the ration of f electrons for the two spin orbit peaks (4f5/2 and 4f7/2) will be
	     	* 6 electrons for 4f5/2 (also denoted N6) and 8 electrons for  4f7/2 (denoted N7). They have different binding energies which can be found in different databases
	     	* The energies are given for the elements in their natural forms relative to the Fermi level for metals
	     	*http://henke.lbl.gov/optical_constants/pert_form.html
	     	*https://srdata.nist.gov/xps/EngElmSrchQuery.aspx?EType=PE&CSOpt=Retri_ex_dat&Elm=Au
	     	*/
	    	double E_4f_min;

	    	//! Start time for odes
	    	double t_start;
	    	//! End time for odes
	    	double t_end;
	    	//! Step size for odes
	    	double t_step;

		//lifetime of 5p-states
		double lifetime_5p_states;
		double rate_5p_states;
	
		//Alias
		double& n_6sp = x[0], & n_5d = x[1], & n_5p = x[2], & n_4f = x[3], & n_val = x[4];
		double& u_val = x[5], & T_e = x[6], & T_ph = x[7];
		double& mu_6sp = x[8], & mu_5d = x[9], & mu_eq = x[10];

};

#endif /* FOURBANDWITHPHONONS_H */


