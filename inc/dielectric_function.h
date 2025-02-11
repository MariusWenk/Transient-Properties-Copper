/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   dielectric_function.h
 * Author: ndione
 *
 * Created on February 9, 2017, 11:17 AM
 */

#ifndef DIELECTRIC_FUNCTION_H
#define DIELECTRIC_FUNCTION_H

#include <gsl/gsl_integration.h>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <stdlib.h> 
#include <cstdlib>
#include <complex>
#include <memory>
#include <utility>

#include "physical_constants.h"
#include "integration.h"
#include "gsl_function_pp.h"
#include "atom_unit_physics.h"
#include "used_data.h"
//#include "fermi_distribution.h"



class dielectric_function {                   

public:

	dielectric_function();
    	virtual ~dielectric_function();
    
    	virtual double GetRealPart(double omega, double Te, double Tph);
	virtual double GetRealPart(double omega, double Te);
	virtual double GetRealPart(double omega, double Te, double Tph, double wavevector);

	virtual double GetImaginaryPart(double omega, double Te, double Tph);
	virtual double GetImaginaryPart(double omega, double Te);
	virtual double GetImaginaryPart(double omega, double Te, double Tph, double wavevector);

	virtual double GetRealFromImaginaryPart(double omega, double Te, double Tph);
	virtual double GetRealFromImaginaryPart(double omega, double Te);
	virtual double GetRealFromImaginaryPart(double omega, double Te, double Tph, double wavevector);
	

	/*************************************************************************************************************************
	* All functions given below need the dielectric function calculated in derived classes
	*Reflectivities are calculated assuming a Three-layer structure vacuum-gold-vacuum
	*If a substrate is used, just replace the refractive index of the gold-substrate interface n_bottom(epsilon) by it's value
	*The three-layer structure can be extended for multiple films. Please use omega in eV
	**************************************************************************************************************************/
 
	//!Optical properties of thin gold film. For bulk material, the same function can be used with larger film thickness
	double NormalIncidenceReflectivityFilm(double omega, std::complex<double> epsilon);
	double sPolarizedReflectivityFilm(double omega, std::complex<double> epsilon); 
	double pPolarizedReflectivityFilm(double omega, std::complex<double> epsilon);
	double NormalIncidenceReflectivityFilm(double omega, double Te, double Tph, std::complex<double> epsilon);
	double sPolarizedReflectivityFilm(double omega, double Te, double Tph, std::complex<double> epsilon); 
	double pPolarizedReflectivityFilm(double omega, double Te, double Tph, std::complex<double> epsilon);

	double NormalIncidenceTransmissivityFilm(double omega, std::complex<double> epsilon);
	double sPolarizedTransmissivityFilm(double omega, std::complex<double> epsilon); 
	double pPolarizedTransmissivityFilm(double omega, std::complex<double> epsilon); 
	double NormalIncidenceTransmissivityFilm(double omega, double Te, double Tph, std::complex<double> epsilon);
	double sPolarizedTransmissivityFilm(double omega, double Te, double Tph, std::complex<double> epsilon); 
	double pPolarizedTransmissivityFilm(double omega, double Te, double Tph, std::complex<double> epsilon); 
 
	double NormalIncidenceAbsorptionFilm(double omega, std::complex<double> epsilon);
	double sPolarizedAbsorptionFilm(double omega, std::complex<double> epsilon); 
	double pPolarizedAbsorptionFilm(double omega, std::complex<double> epsilon);
	double NormalIncidenceAbsorptionFilm(double omega, double Te, double Tph, std::complex<double> epsilon);
	double sPolarizedAbsorptionFilm(double omega, double Te, double Tph, std::complex<double> epsilon); 
	double pPolarizedAbsorptionFilm(double omega, double Te, double Tph, std::complex<double> epsilon);

	// Fields reflectivity and transmissivity --> Fresnel coeffs they are necessary to calculate the optics
	// Top for upper layer, bottom for down layer. 
	std::complex<double> r_top_film_normalIncidence(std::complex<double> epsilon); // normal incidence
	std::complex<double> r_top_film_sPolarized(std::complex<double> epsilon); //s-pol
	std::complex<double> r_top_film_pPolarized(std::complex<double> epsilon); //p-pol
	std::complex<double> t_top_film_normalIncidence(std::complex<double> epsilon); // normal incidence
	std::complex<double> t_top_film_sPolarized(std::complex<double> epsilon); //s-pol
	std::complex<double> t_top_film_pPolarized(std::complex<double> epsilon); //p-pol

	std::complex<double> r_film_bottom_normalIncidence(std::complex<double> epsilon); // normal incidence
	std::complex<double> r_film_bottom_sPolarized(std::complex<double> epsilon); //s-pol
	std::complex<double> r_film_bottom_pPolarized(std::complex<double> epsilon); //p-pol
	std::complex<double> t_film_bottom_normalIncidence(std::complex<double> epsilon); // normal incidence
	std::complex<double> t_film_bottom_sPolarized(std::complex<double> epsilon); //s-pol
	std::complex<double> t_film_bottom_pPolarized(std::complex<double> epsilon); //p-pol

	//Refractive index considered medium: Three-layer structure
	std::complex<double> n_top(std::complex<double> epsilon);
	std::complex<double> n_bottom(std::complex<double> epsilon);
	std::complex<double> n_film(std::complex<double> epsilon);


	double TransmittedAngle();
	std::complex<double>cosineRefractedAngle(std::complex<double> epsilon);
	//Function which take into account the film Thickness
	std::complex<double> BetaFunction(double omega, std::complex<double> epsilon);
	//Get private variables
	double GetAngleIncidence();
	double GetFilmThickness();
	double GetDielectricConst();
	
	//Single photon absorption
	double SinglePhotonAbsorption(double omega, std::complex<double> epsilon);
	//Material electrical conductivity
	std::complex<double> ComplexConductivity(double omega, std::complex<double> epsilon);
	double RealConductivity(double omega, std::complex<double> epsilon);
	double bulkReflectivity(std::complex<double> epsilon);
	double ImaginaryConductivity(double omega , std::complex<double> epsilon);
	double SkinDepth(double omega, std::complex<double> epsilon);
 
	//General case solved with transfer matrix method();
	std::complex<double>MatrixElement_m11(std::complex<double> epsilon, bool polarization);
	std::complex<double>MatrixElement_m12(std::complex<double> epsilon, bool polarization);
	std::complex<double>MatrixElement_m21(std::complex<double> epsilon, bool polarization);
	std::complex<double>MatrixElement_m22(std::complex<double> epsilon, bool polarization);
	
	double normalIncidenceReflectivityBulk(std::complex<double>);

protected:

	double angle_incidence ;  
	double film_thickness; 
	double dielectric_const_value;  
	Cintegration* integrator; // Integration Class pointer
    

};

#endif /* DIELECTRIC_FUNCTION_H */