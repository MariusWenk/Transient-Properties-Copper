/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   D_band.cpp
 * Author: ndione
 * 
 * Created on 4. April 2018, 11:41
 */

#include <memory>
#include <random>

#include "D_band.h"

D_band::D_band() {
	
	this->d_band_mass = d_band_mass_const*m_e_au;  
    	this ->unit_cell_volume = unit_volume_au;
	this ->n_d_density = number_of_d_electrons / unit_cell_volume;
	this ->fermi_energy = fermi_energy_const *eV_to_Hartree;
	this ->distance_d_band_E_F = distance_d_band_E_F_const *eV_to_Hartree;	 
}

D_band::D_band(int number_disc_points) {

    	this ->number_disc_points = number_disc_points;
}

D_band::~D_band() {

}

double D_band::Get_d_BandMass() {

    	return d_band_mass;
}

int D_band::GetNumberDiscPoints() {

    	return number_disc_points;
}

double D_band::Get_n_d_Density() {

    	return n_d_density;
}

double D_band::GetUnitVolume() {

    	return unit_cell_volume;
}

double D_band::GetFermiEnergy() {

	return fermi_energy;
}

double D_band::GetDistance_d_BandFermi() {

	return distance_d_band_E_F;
}

double D_band::Get_d_DOS(double energy) {
 
	if ((energy >= 0) && (energy < GetFermiEnergy() - GetDistance_d_BandFermi())) {

	        return (std::sqrt(2.0 * (energy - 0.0)) * std::pow(Get_d_BandMass(), 1.5)) / (PI * PI * std::pow(hbar_au, 3.0));
	}

	return 0.0;     
}


double D_band::GetDispersionRelation(double wavevector) {

    	return hbar_au * hbar_au * wavevector * wavevector / (2.0 * Get_d_BandMass());
}

double D_band::Get_d_DOSWavevector(double wavevector) {

    	return Get_d_DOS(GetDispersionRelation(wavevector));
}

double D_band::GetWavevector(double energy) {

    	return std::sqrt(2.0 * Get_d_BandMass() * GetDispersionRelation(energy) * GetDispersionRelation(energy) / hbar_au);

}
