/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   D_band.h
 * Author: ndione
 *
 * Created on 4. April 2018, 11:41
 */



#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h> 


#include "physical_constants.h"
#include "atom_unit_physics.h"
#include "BoostRootFinding1D.h"

#ifndef D_BAND_H
#define D_BAND_H

class D_band {

public:

    	// Constructor
    	D_band();
    	// Copy constructor
    	D_band(int number_disc_points);
    	// Destructor
    	virtual ~D_band();

	//! Hands over the DOS at a given energy
    	double Get_d_DOS(double energy);
    	//! Hands over the DOS at a given number of discretization points
    	double Get_d_DOS(int number_of_disc_points);
   	//! Hands over the DOS at a given wavevector 
    	double Get_d_DOSWavevector(double wavevector);
    	//! Hands over the band mass
    	double Get_d_BandMass();
    	//! Hands over the unit cell volume of the material
    	double GetUnitVolume();
    	//! Hands over electron density of d band
    	double Get_n_d_Density();
    	//! hands over the number of discretization points
    	int GetNumberDiscPoints();
    	//! Hands over fermi energy of the material
    	double GetFermiEnergy();
	//!Hands over distance from the d-band edge to the Fermi Energy
	double GetDistance_d_BandFermi();
    	//! Hands over the dispersion relation 
    	double GetDispersionRelation(double wavevector);
    	//! Hands over the wavevector
    	double GetWavevector(double energy);

private:

	// Effective mass of sp band
    	double d_band_mass;
    	// Number of discretization points
    	int number_disc_points;
    	// Unit cell volume
    	double unit_cell_volume;
    	// Free electron density
    	double n_d_density;
    	// Fermi energy
    	double fermi_energy;
   	//distance d-band-edge to Fermi Energy
	double distance_d_band_E_F; 
   

};

#endif /* D_BAND_H */

