/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Laser.h
 * Author: Ndione, Ag Rethfeld
 *
 * Created on January 24, 2018, 10:44 AM
 */

#ifndef LASER_H
#define LASER_H

#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <stdlib.h> 
#include <complex>

#include "physical_constants.h"
#include "atom_unit_physics.h"
#include "used_data.h"
#include "dielectric_function.h"
#include "dielectric_function_drude.h"
#include "dielectric_function_lorentz.h"
#include "BoostRootFinding1D.h"


/*!********************************************************
 * \brief Contains properies of the laser pulse
 *
 * Contains pulse duration, fluence and photon energy of
 * the laser pulses. Calculates intensity directly beneath 
 * the surface.
 **********************************************************/

class Claser {
   

public:
	//! Constructor
	Claser();
	Claser(bool boolean);
    	//! Destructor
    	~Claser();


    	//! Hands over pulse duration
    	double GetPulseDuration();
    	//! Hands over fluence  
    	double GetFluence();
    	//! Hands over photon energy  
    	double GetPhotonEnergy();
    	//! Hands over wavelength 
    	double GetWavelength();
    	//! Hands over angular frequency of the laser pulse
    	double GetOmega();
    	//! Hands over electric field 
    	double LaserElectricField_eq(double time, double Te, double Tph);
    	double LaserElectricField_neq(double time, double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d);
   	 //! Hands over the magnetic field 
    	double LaserMagneticField_eq(double time, double Te, double Tph);
    	double LaserMagneticField_neq(double time, double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d);
    	//Hands over a Gaussian source term
    	double LaserSourceTerm(double time);
	//! Hands over Gaussian laser intensity
    	double LaserIntensity(double time, double Te, double Tph);
	double LaserIntensity_neq(double time, double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d);
	//! Hands over the time at which the Gaussaian pulse has its maximum
	double GetTimePeak();
	//!Hands over the reflectivity	
	double CalcReflectivity(double Te, double Tph);
	double CalcReflectivity_neq(double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d);
	double CalcTransmissivity(double Te, double Tph);
	double CalcTransmissivity_neq(double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d);
	double CalcAbsorption(double Te, double Tph);
	double CalcAbsorption_neq(double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d);
	//!Hands over the total absorbed energy in J/kg
	double GetAbsorbedEnergy();
	//!Hands over the material density
	double GetMaterialDensity();
	//! hands over the absorption cefficient alpha
	double CalcAbsorptionCoeff(double Te, double Tph);
	double CalcAbsorptionCoeff_neq(double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d);
	//!Hands over the  total dielectric function
	std::complex<double> CalcDielectricFunction(double Te, double Tph);
	std::complex<double> CalcDielectricFunction_neq(double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d);	
	//!Hands over the Keldysh parameter
	double KeldyshParameter(double);
	double KeldyshParameter();
	//!Hands over incident source term
	double IncidentSourceTerm(double time, double Te, double Tph);
	double IncidentSourceTerm_neq(double time, double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d);
	//!Hands over Incident laser intensity
	double IncidentLaserIntensity(double time, double Te, double Tph);
	double IncidentLaserIntensity_neq(double time, double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d);

	
  
private:

    	//! Time for peak Gaussian pulse;
	double time_peak; 
    	//! Fluence in J/m<sup>2</sup>
    	double fluence;
    	//! Pulse duration  
    	double pulse_duration;
    	//! Laser wavelength  
    	double wavelength;
    	//! Photon energy  
    	double photon_energy;
    	//! Angular frequency of the laser pulse 
    	double omega;
    	//!Maximum Bessel order that should be used
    	int max_bessel_order;
 
    	double absorbed_energy;
	double material_density;

	dielectric_function *diel_pointer;
	dielectric_function_drude *drude_pointer;
	dielectric_function_lorentz *lorentz_pointer;
	BoostRootFinding1D *boost_pointer;
	bool boolean;
 

	  
};



#endif /* LASER_H */

