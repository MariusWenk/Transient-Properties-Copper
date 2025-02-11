/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   dielectric_function_drude.h
 * Author: ndione, Ag Rethfeld
 *
 * Created on February 9, 2017, 11:18 AM
 */

#ifndef DIELECTRIC_FUNCTION_DRUDE_H
#define DIELECTRIC_FUNCTION_DRUDE_H

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
#include "BoostRootFinding1D.h"



class dielectric_function_drude : public dielectric_function {
    
public:
	dielectric_function_drude();
	dielectric_function_drude(BoostRootFinding1D *boost_point);
    	virtual ~dielectric_function_drude();

       
 	//Complex Drude dielectric function
    	std::complex<double> ComplexDielectricDrude(double omega, double Te, double Tph);
	std::complex<double> ComplexDielectricDrude(double omega, double Te, double Tph, double wavevector);
	//Complex Drude dielectric function for noneq band occupation 
	std::complex<double> ComplexDielectricDrude_neq(double omega, double Tph, double sp_density, double d_density);
	// Hands over the real part 
	double GetRealPart(double omega, double Te, double Tph);
	double GetRealPart(double omega, double Te, double Tph, double wavevector);
	double GetRealPart_neq(double omega, double Tph, double sp_density, double d_density);
	// Hands over the imaginary part
	double GetImaginaryPart(double omega, double Te, double Tph);
	double GetImaginaryPart(double omega, double Te, double Tph, double wavevector);
	double GetImaginaryPart_neq(double omega, double Tph, double sp_density, double d_density);
	//Hands over the imaginary part using Kramers-Kronig
	double GetRealFromImaginaryPart(double omega, double Te, double Tph);
	double GetRealFromImaginaryPart(double omega, double Te, double Tph, double wavevector);
	double GetRealFromImaginaryPart_neq(double omega, double Tph, double sp_density, double d_density);
    	    	 
    	 
    	//e-ph collision frequency
	double ElectronPhononCollisionFreq(double ph_temp); 
	//e-e collision frequency  	
	double ElectronElectronCollisionFreq(double Te);
	//e-e collision frequency  for noneq band occupation
	double ElectronElectronCollisionFreq_neq(double d_density);
	//total collision frequency
	double TotalCollisionFreq(double Te, double Tph);
	//total collision frequency for noneq. band occupation
	double TotalCollisionFreq_neq(double Tph, double d_density);
	
	//!Hands over Fermi velocity
	double FermiVelocity();   
	// Hands over Te-dependent dc-conductivity
	double dcConductivity(double Te, double Tph); 
	double dcConductivity_neq(double Tph, double sp_density, double d_density);
	// Hands over Te-dependent resiytivity
	double Resistivity(double Te, double Tph);
	double Resistivity_neq(double Tph, double sp_density, double d_density);
	// Collision time	
	double GetTotalCollisionTime(double Te, double Tph);
	//plasma freq	
	double GetPlasmaFreqHighTe(double Te); 
	double CalcPlasmaFreqHighTe(double sp_density); 
	double GetPlasmaFreqHighTeSquare(double Te);
	double CalcPlasmaFreqHighTeSquare(double sp_density);

	//Hands over private variables
	double GetDielectricConst();
    	double GetTemperatureStart();   
    	double Get_n_sp_eq();
    	double Get_n_d_eq();   	
    	double Get_N_sp_Equilibrium();
    	double Get_N_d_Equilibrium();
   	double GetFermiEnergy();
    	double GetUnitVolume();
    	double GetBandMass_sp();
	double GetBandMass_d();
	double get_el_el_const();
	double get_el_ph_coll_freq_cold();
	
    

private:

    	double unit_volume;
    	double dielectric_const_value;  	
    	double temperature_start;
    	double n_sp_eq;
    	double n_d_eq;
    	double N_sp_eq;
    	double N_d_eq;
    	double fermi_energy; 
    	double sp_band_mass;
    	double d_band_mass;

     	double el_el_const;
     	double el_ph_coll_freq_cold;
     	
    	BoostRootFinding1D  *boost_point; 
    	Cintegration *int_pointer; 

};

#endif /* DIELECTRIC_FUNCTION_DRUDE_H */