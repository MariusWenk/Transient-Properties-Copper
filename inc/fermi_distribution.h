 
/* 
 * File:   fermi_distribution.h
 * Author: Ndione, Ag Rethfeld
 *
 * Created on 22. Februar 2018, 09:48
 */

#ifndef FERMI_DISTRIBUTION_H
#define FERMI_DISTRIBUTION_H

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <functional>
#include <fstream>
#include <array>
#include <string>
#include <complex>
#include <limits>
#include <tuple>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_const.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_const.h>
#include <boost/math/tools/roots.hpp>


#include "atom_unit_physics.h"
#include "physical_constants.h"
#include "used_data.h"
#include "gsl_function_pp.h"
#include "integration.h"
#include "root_finding.h"
#include "sp_band.h"
#include "D_band.h"
#include "dielectric_function.h"
#include "dielectric_function_drude.h"

class Cfermi_distribution {

public:

    
    	std::vector<std::vector<double>> dosfile; // Matrix or vector of vectors which will contain all 4 columns (E, sp_dos, d_dos, sp+d_dos)
    	//Array to store Te and mu_eq
    	std::vector<std::vector<double>> test_file;

    	//! Constrcutor
    	Cfermi_distribution();
    	//! Copy constructor
    	Cfermi_distribution(std::string dos_filename);
    	//! Destrcutor
    	virtual ~Cfermi_distribution();


    	//Read a file and hands over energy needed density
    	void ReadDOSFiles(const char *filename, std::vector<double> &E, std::vector<double> &SP_dos, std::vector<double> &D_dos, std::vector<double> &SPD_dos);
    	//! Hands over the fermi energy
    	double GetFermiEnergy();
	//!Hands over starting electron temperature
	double GetTemperatureStart();
    	//! Hands over the number of dicretization points
    	double GetNumberDiscPoints();
    	//! Hands over the number of states available for absorption in the sp-band in the whole DOS range. This accounts also Pauli blocking
    	double NumberStatesTotal_sp(double T, double mu_sp, double mu_d, double photon_energy);
    	//! Hands over the number of states available for absorption in the d-band in the whole DOS range. This accounts also Pauli blocking
    	double NumberStatesTotal_d(double T, double mu_sp, double mu_d, double photon_energy);
    	double NumberStatesTotal(double T, double mu_sp, double mu_d, double photon_energy);
    	//! Hands over interband percentage  in the whole DOS range
    	double InterbandPercentageTotal(double T, double mu_sp, double mu_d, double photon_energy);
    	//! Hands over intraband percentage  in the whole DOS range
    	double IntrabandPercentageTotal(double T, double mu_sp, double mu_d, double photon_energy);
    	//! Hands over the fermi distribution
	double GetFermiDistribution(double energy, double T, double mu);
    	//! Hands over the internal energy
    	double CalcInternalEnergy(double T, double mu);
    	//! Hands over the internal energy of sp electrons
    	double CalcInternalEnergy_sp(double T, double mu);
    	//! Hands over the internal energy of d electrons
    	double CalcInternalEnergy_d(double T, double mu);
    	//! Hands over the electrons density
    	double CalcDensity(double T, double mu);
    	//! Hands over equilibrium density for sp electrons: this function is also needed for the three and two-band models
    	double CalcEquilibriumDensity_sp(double T, double mu);
    	//! Hands over equilibrium for d electrons
    	double CalcEquilibriumDensity_d(double T, double mu);
    	//! Hands over the chemical_potential using gsl root finding
    	double GetChemicalPotential_GSL(double temperature);
     	//! Hands over chemical potential using Boost root finding
    	double CalcMuRootFindingBoost(double n0, double T, double E_Fermi);
    	//! Hands over the partial heat capacity CT=\partial u/\partial Te
    	//! Do not confuse with total heat capaciaty Ce = du/dT
    	double Get_CT(double T, double mu);
    	//! Hands over the partial heat capacity  of sp electrons CTsp=\partial usp/\partial Te
    	double Get_CT_sp(double T, double mu);
    	//! Hands over the partial heat capacity  of d electrons CTd=\partial ud/\partial Te
    	double Get_CT_d(double T, double mu);
    	//! Hands over the partial heat capacity  Cmu = \partial u/\partial\mu
    	double Get_CT_mu(double T, double mu);
   	//! Hands over the partial heat capacity  Cmu_sp = \partial usp/\partial\musp
    	double Get_CT_mu_sp(double T, double mu);
    	//! Hands over the partial heat capacity  Cmu_d = \partial ud/\partial/mud
    	double Get_CT_mu_d(double T, double mu);
    	//! Hands over the variation of the density with respect to temperature dn_dTe
    	double Get_P(double T, double mu);
    	//! Hands over the variation of sp density with respect to temperature dn_sp_dTe
    	double Get_P_sp(double T, double mu);
    	//! Hands over the variation of d density with respect to temperature dn_d_dTe
    	double Get_P_d(double T, double mu);
    	//! Hands over the variation of the density with respect to chemical_potential dn_dmu
    	double Get_P_mu(double T, double mu);
    	//! Hands over the variation of sp density with respect to chemical_potential dnsp_dmu
    	double Get_P_mu_sp(double T, double mu);
    	//! Hands over the variation of d density with respect to chemical_potential dnd_dmu
    	double Get_P_mu_d(double T, double mu);
    	//! total heat capacity
    	double heatCapacity(double, double);
    	//!sp-heat capacity
    	double heatCapacity_sp(double, double);
    	//!d-heat capacity
    	double heatCapacity_d(double, double);
    	//! Hands over the transient entropy of the system
    	double GetEntropy(double time, double energy);
    	//! Hand over relaxation time according to Fermi-liquid theory
    	double RelaxationTime(double Te, double energy);
    	//! Hands over total absorption coefficient using the total absorbed energy from the laser
    	double AlphaTotFromEnergyDensity(double time);
    	//! Lower bond used in the GSL integration routine
    	double LowerBound();
    	//!Upper bond for the GSL integration routine
    	double UpperBound();    
    	//! Read The dosfile, makes an interpolation and return sp-band DOS
    	double GetDOS_sp(double energy);
    	//! Read The dosfile, makes an interpolation and return d-band DOS
    	double GetDOS_d(double energy);
    	//! Read The dosfile, makes an interpolation and return total DOS
    	double GetDOSTot(double energy);
    	//! Hands over the derivative of the Fermi distribution with respect to electrons temperature
    	double df_dTe(double energy, double T, double mu);
    	//! Hands over the derivative of the Fermi distribution with respect to the chemical potential
    	double df_dmu(double energy, double T, double mu);
	double GetDistance_d_BandtoFermi();
	double StepFunction(double energy);
	double SmoothFermiFunction(double energy);
	double GetUnitVolume();
	double TestGSL();

	 
private:

    	//GSl accelerator of sp-dos
    	gsl_interp_accel *sp_dos_accel;
    	//GSL accelerator of d-dos
    	gsl_interp_accel *d_dos_accel;
    	//GSL accelerator of tot-dos
    	gsl_interp_accel *tot_dos_accel;
    	//GSL spline of sp-electrons
    	gsl_spline *sp_dos_spline;
    	//GSL spline of d-electrons
    	gsl_spline *d_dos_spline;
    	//GSL spline of tot-electrons
    	gsl_spline *tot_dos_spline; 	

    	// Array for energy
    	double* energy;
    	// Array for total DOS
    	double* DOS;
    	// Array for sp DOS
    	double *DOS_sp;
    	//Array for d DOS
    	double *DOS_d;
 

    	// Fermi energy
    	double fermi_energy;
    	// Number of discretization points
    	unsigned int n;
    	// Class pointer to sp_band Class
    	Sp_band *sp_band_pointer;
    	D_band *d_band_pointer;
    	unsigned int dos_factor;
    	double *laser_ptr;

    	Cintegration *integrator_ptr;
    	 
	//! distance from d-band edge to Fermi-energy	
	double distance_d_band_E_F;
 	double temperature_start;

	double unit_volume;

};


struct root_params {
     
};


#endif /* FERMI_DISTRIBUTION_H */

