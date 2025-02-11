/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   dielectric_function_lindhard.cpp
 * Author: ndione, Ag Rethfeld
 * 
 * Created on February 27, 2017, 4:18 PM
 */


#include"dielectric_function_lindhard.h"
#include "dielectric_function_lindhard.h"

/*
dielectric_function_lindhard::dielectric_function_lindhard() {
  
    integrator = new Cintegration;
}

dielectric_function_lindhard::~dielectric_function_lindhard() {
    	
	delete integrator;
}

void dielectric_function_lindhard::do_timestep() {

    double prefactor = pow(electron_charge, 2.0) / vacuum_permittivity / 4 / M_PI / M_PI;

    auto integrand = [ = ] (double energy) -> double {

        return band->get_DOS(energy) * band->get_distribution(energy) / band->get_bandmass(energy);
    };

    integral_optical_limit = prefactor * integrator->integrate_qag(integrand, band->get_E_min(), band->get_E_max());
}

// Calculate the Real part of Lindhard function in the optical limit (q-->0)

//double dielectric_function_lindhard::CalcRealPartEquilibrium(double omega) {
//    return 1.0 - integral_optical_limit / omega / omega;
//}

//double dielectric_function_lindhard::CalcRealPartEquilibrium(double omega, double wavenumber_q) {
//    if (wavenumber_q == 0.0) {
//        return CalcRealPartEquilibrium(omega);
//    }
//    else {
//        exit(1);
//    }
//}

// Calculate Imaginary part of Lindhard in the optical limit with Kramers-Kronig 

double dielectric_function_lindhard::CalcImaginaryPartEquilibrium(double omega) {

    auto integrand = [ = ] (double x) -> double {

        return x * x * (CalcRealPartEquilibrium(x) - 1) / (x + omega);
    };

    return -(2.0 / PI / omega) * integrator -> integrate_qawc(integrand, 10., 1.0e25, omega);
}
double dielectric_function_lindhard::CalcImaginaryPartEquilibrium(double omega, double wavevector){

    return -1;
}
double dielectric_function_lindhard::CalcRealPartEquilibrium(double omega){

    return -1;
}

double dielectric_function_lindhard::CalcRealPartEquilibrium(double omega, double wavenumber_q) {
    
if (wavenumber_q == 0) {

        return CalcImaginaryPartEquilibrium(omega);

    } else {

        double result = (4.0 * pow(PI, 4.0) * electron_charge * electron_charge) / (band->get_spin_factor() * pow(wavenumber_q, 3.0));

        auto integrand = [ = ] (double energy) -> double {

            double energy_prime = energy + hbar*omega;

            return energy_prime*result;// this is not correctly implemented yet
        };
		return integrator->integrate_qag(integrand, 0,10);// This is not correctly implemented yet
    }
}
*/


//void dielectric_function_lindhard::do_timestep() {
//    
//    double prefactor = pow(electron_charge, 2.0) / vacuum_permitivity;
//    auto integrand = [ = ] (double energy) -> double {
//        return band->get_DOS(energy) * band->get_distribution(energy) / band->get_bandmass(energy);
//    };
//    integral_optical_limit = prefactor * integrator->integrate_qag(integrand, band->get_E_min(), band->get_E_max());
//}
