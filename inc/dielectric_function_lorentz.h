/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   dielectric_function_lorentz.h
 * Author: ndione, Ag Rethfeld
 *
 * Created on 27. März 2018, 15:13
 */

#ifndef DIELECTRIC_FUNCTION_LORENTZ_H
#define DIELECTRIC_FUNCTION_LORENTZ_H


#include <gsl/gsl_integration.h>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <stdlib.h> 
#include <cstdlib>
#include <complex>

#include"dielectric_function.h"
#include "physical_constants.h"
#include "integration.h"
#include "gsl_function_pp.h"
#include "atom_unit_physics.h"
#include "fermi_distribution.h"
#include "BoostRootFinding1D.h"
#include "used_data.h"



class dielectric_function_lorentz : public dielectric_function {
    
public:
    	//! Constructor
	dielectric_function_lorentz();
	dielectric_function_lorentz(BoostRootFinding1D  *boost_point);
    	//! Destructor 
  	virtual ~dielectric_function_lorentz();

    	//! Hands over complex interband dielactric function
    	std::complex<double> ComplexDielectricLorentz(double omega, double Te);
	std::complex<double> ComplexDielectricLorentz_neq(double omega, double Te, double mu_sp, double mu_d, double d_density);
	// Real part of Lorentz dilectric function
    	double GetRealPart(double omega, double Te);
	double GetRealPart_neq(double omega, double Te, double mu_sp, double mu_d, double d_density);
	//Imaginary part of Lorentz dielectric function
	double GetImaginaryPart(double omega, double Te);
	double GetImaginaryPart_neq(double omega, double Te, double mu_sp, double mu_d, double d_density);
	//Imaginary part that will be calculated with Kramers-kronig 
	double GetRealFromImaginaryPart(double omega, double Te);
	double GetRealFromImaginaryPart_neq(double omega, double Te, double mu_sp, double mu_d, double d_density);

    
    	//! Hands over Lorentz oscillators strength at room Te
    	double GetEqStrength1();
    	double GetEqStrength2();
    	double GetEqStrength3();
    	double GetEqStrength4();
    	double GetEqStrength5();
	//! Hands over Lorentz oscillators strength at elevated Te
	double GetNonEqStrength1(double Te);
	double GetNonEqStrength2(double Te);
	double GetNonEqStrength3(double Te);
	double GetNonEqStrength4(double Te);
	double GetNonEqStrength5(double Te);
	//! Hands over Lorentz oscillators strength at elevated Te for noneq. band occupation
	double GetNonEqStrength1_neq(double, double);
	double GetNonEqStrength2_neq(double, double);
	double GetNonEqStrength3_neq(double, double);
	double GetNonEqStrength4_neq(double, double);
	double GetNonEqStrength5_neq(double, double);

    	//! Hands over interband transition frequencies at room Te
    	double GetOmega1();
    	double GetOmega2();
    	double GetOmega3();
    	double GetOmega4();
    	double GetOmega5();

	//! Hands over interband transition frequencies at elevated Te
	double OmegaNonEq1(double);
	double OmegaNonEq2(double);
	double OmegaNonEq3(double);
	double OmegaNonEq4(double);
	double OmegaNonEq5(double);
	// case of nonequilibrium band occupation
	double OmegaNonEq1_neq(double, double);
	double OmegaNonEq2_neq(double, double);
	double OmegaNonEq3_neq(double, double);
	double OmegaNonEq4_neq(double, double);
	double OmegaNonEq5_neq(double, double);

	//!Hands over interband collision frequencies at room Te
    	double GetGamma1();
    	double GetGamma2();
    	double GetGamma3();
    	double GetGamma4();
    	double GetGamma5();

	//!Hands over interband collision frequencies at elevated Te
	double GammaNonEq1(double Te);
	double GammaNonEq2(double Te);
	double GammaNonEq3(double Te);
	double GammaNonEq4(double Te);
	double GammaNonEq5(double Te);
	//!Hands over interband collision frequencies at elevated Te in case of noneq band occupation
	double GammaNonEq1_neq(double d_density);
	double GammaNonEq2_neq(double d_density);
	double GammaNonEq3_neq(double d_density);
	double GammaNonEq4_neq(double d_density);
	double GammaNonEq5_neq(double d_density);

    	// particle density and particle number at room-Te
    	double Get_n_sp_eq();
    	double Get_n_d_eq();  	
	double Get_N_sp_eq();
	double Get_N_d_eq();
	
	//Hands over the Fermi energy
	double GetFermiEnergy();
	//Hands over the starting temperature
	double Get_Te_startTemperature();
	//Hands over the volume of the unit cell... Au here
	double GetUnitVolume();
	
	//Hands over complex Lorentz functions at elevated Te
 	std::complex<double> ComplexLorentzOscillator_1(double omega, double Te);
	std::complex<double> ComplexLorentzOscillator_2(double omega, double Te);
	std::complex<double> ComplexLorentzOscillator_3(double omega, double Te);
	std::complex<double> ComplexLorentzOscillator_4(double omega, double Te);
	std::complex<double> ComplexLorentzOscillator_5(double omega, double Te);
	//Hands over complex Lorentz functions for nonequilibirum band occupation
	std::complex<double> ComplexLorentzOscillator_1_neq(double omega, double Te, double mu_sp, double mu_d, double d_density);
	std::complex<double> ComplexLorentzOscillator_2_neq(double omega, double Te, double mu_sp, double mu_d, double d_density);
	std::complex<double> ComplexLorentzOscillator_3_neq(double omega, double Te, double mu_sp, double mu_d, double d_density);
	std::complex<double> ComplexLorentzOscillator_4_neq(double omega, double Te, double mu_sp, double mu_d, double d_density);
	std::complex<double> ComplexLorentzOscillator_5_neq(double omega, double Te, double mu_sp, double mu_d, double d_density);


	//!Hands over functions used for the amplitudes of the Lorentz oscillators at elevated Te
	double relativeDistribution(double energy, double Te);
	//possible approximation for Lorentz damping at elevated Te
	double dampingApproximation(double);
	double dampingApproximation_neq(double); 
	//possible approximation for Lorentz interband frequencies at elevated Te
	double OmegaApproximation(double);
	double OmegaApproximation_neq(double, double);
	
	//return const for electron-electron scattering rate
	double get_const_el_el();



private:

    	double unit_volume;
	double fermi_energy;
	double Te_start_temperature;
	double n_sp_eq;
    	double n_d_eq;
	double N_sp_eq;
	double N_d_eq;
	int oscillator_number;

    	// Oscillator stenght at equilibrium: 300 K
		std::array<double,osci_num_const> eq_strength;

    	// Interband transitions frequencies
		std::array<double,osci_num_const> omega;

    	// 1/Collision time for interband transition at room temperature
		std::array<double,osci_num_const> gamma;
	
	//const for electron-electron scattering rate
	double const_el_el;

  	Cintegration *int_pointer;
  	BoostRootFinding1D  *boost_point;
};

#endif /* DIELECTRIC_FUNCTION_LORENTZ_H */
