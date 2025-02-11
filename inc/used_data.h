/*
 * File:   used_data.h
 * Author: Wenk, Ag Rethfeld
 *
 * Created on 16. May 2022, 11:21
 */

#include <math.h>

#ifndef USED_DATA_H
#define USED_DATA_H

// The constants defined here are only prefactors and the transformations to the right units are given when called. This is because different units may be used.

// Free parameters:     BoostRootFinding1D::GetMuRootFinding(double T) -> Newton Raphson Borders
// (not defined here    Some values in unused (depricated) constructors, classes or fucntions
// but elsewhere)       Absorption coefficients: alpha_tot

//                      dielectric_function_lorentz::dielectric_function_lorentz(BoostRootFinding1D  *boost_point) : dielectric_function() -> omega, gamma, strength_eq values
//                      transient_optical_param::transient_optical_param(const char *filename,  BoostRootFinding1D *boost_point) -> omega, gamma, strength_eq values

/*************************************************/
// Data for Gold:

//#define DOS_input ("input/DOS_Au/SP_D_SPD_states_EH_atom.in")

//#define e_ph_coupling ("input/coupling_au/Lin_coupling_au.in")
//#define e_ph_coupling ("input/coupling_au/Held_coupling_eq_au.in")
//#define e_ph_coupling ("input/coupling_au/Smirnov_coupling_au.in")
//#define e_ph_coupling ("input/coupling_au/Petrov_coupling_au.in")
//#define e_ph_coupling ("input/coupling_au/Medvedev_coupling_au.in")

#define e_ph_coupling_sp ("input/coupling_au/Held_sp_coupling_au.in") 	//->density resolved sp-coupling only needed for 2band2T model
#define e_ph_coupling_d ("input/coupling_au/Held_d_coupling_au.in") 	//->density resolved  d-coupling only needed for 2band2T model

// #define sp_DOS_stretching (0.9122) //unitless
// #define d_DOS_stretching (1.1231) //unitless
// #define distance_d_band_E_F_const (2.1) //eV
#define sp_band_mass_const (0.6) //number of electron masses (m_e)
#define d_band_mass_const (5.4) //number of electron masses (m_e)
// #define sp_mass_const (0.99) //number of electron masses (m_e)

// #define fermi_energy_const (10.224) //eV

// #define material_density_const (19.3e3) //kg/m³

// #define unit_volume_SI (1.695e-29) //m³ //from Seb DOS

// #define number_of_d_electrons (10.0) //unitless
// #define number_of_sp_electrons (1.0) //unitless

// // Fit parameters Drude-Lorentz:
// #define dielectric_infty (3.104) //unitless
// #define osci_num_const (5)
// // other parameters defined in "dielectric_function_lorentz.cpp"

// #define phonon_heat_capacity_const (2.327e6) //J/(m³*K)

// #define melting_temp_const (1337.33) //K

// #define el_el_collision_const (0.36e15) //s^(-1) //=A/V² //from Fourment et al.
// #define el_ph_collision_const_cold (1.0/11.9e-15) //s^(-1) // for Au from Ng. PRE 2016

//Three Band and Four Band:
#define dist_sp_f_const ((87.6+83.9)*0.5) //eV //83.9 for 4f7/2, 87.6 for 4f5/2
#define dist_sp_p_const ((57.2+74.2)*0.5) //eV
#define lifetime_p_states_const (10e-15) //s
#define d_edge_const (1.9) //eV

#define number_of_p_electrons (6.0) //unitless
#define number_of_f_electrons (14.0) //unitless

#define inherent_linewidth_const ((300e-3+280e-3)*0.5) //eV //280e-3 for 4f_5/2, 300e-3 for 4f_7/2
/*************************************************/



/*************************************************/
// Data for Copper:

//#define DOS_input ("input/DOS_Cu/Cu_l_proj_dos.in")
#define DOS_input ("input/DOS_Cu/DOS_easy_Cu.in")

#define e_ph_coupling ("input/coupling_Cu/Lin_coupling_Cu_au.in") //https://compmat.org/electron-phonon-coupling/

#define sp_DOS_stretching (0.7573) //unitless
#define d_DOS_stretching (1.0372) //unitless
//#define sp_DOS_stretching (1) //unitless
//#define d_DOS_stretching (1) //unitless
#define distance_d_band_E_F_const (2.1) //eV //Obergfell
#define sp_mass_const (1.49) //number of electron masses (m_e) //https://journals.aps.org/prb/pdf/10.1103/PhysRevB.6.4370

#define fermi_energy_const (9.361098) //eV //mu from DOS //difference between smallest non-zero value and Fermi energy at E = 0eV

#define material_density_const (8.94e3) //kg/m³ //https://www.copper.org/resources/properties/atomic_properties.html

//#define unit_volume_SI (1.182e-29) //m³ //https://www.copper.org/resources/properties/atomic_properties.html
#define unit_volume_SI (1.18e-29) //m³ //calculated from atomic mass and density

#define number_of_d_electrons (10.0) //unitless
#define number_of_sp_electrons (1.0) //unitless

#define phonon_heat_capacity_const (3.503e6) //J/(m³*K)
#define el_el_collision_const (0.36e15) //s^(-1) //=A/V² //from Fourment et al.
#define el_ph_collision_const_cold (1.0/6.9e-15) //s^(-1) //https://journals.aps.org/prb/pdf/10.1103/PhysRevB.6.4370
#define melting_temp_const (1358) //K //irrelevant

//Fit parameters Drude-Lorentz:
#define dielectric_infty (3.65223607643857) //unitless
//#define dielectric_infty (3.75624276813676) //unitless
#define osci_num_const (4)
//The parameters are defined in "dielectric_function_lorentz.cpp"
/*************************************************/

/*************************************************/
// Genereal Data:

#define file_name_appendix ("4easy")

#define temperature_start_const (300.0) //K
#define temperature_smooth_fermi_function (1000.0) //K

#define relaxation_time_const (258e-15) //s
//#define relaxation_time_const (400e-15) //s

#define simulation_t_start (-1e-12) //s
//#define simulation_t_start (-500e-15) //s
#define simulation_t_end (10e-12) //s
#define simulation_t_step (1e-16) //s //exact
//#define simulation_t_step (1e-15) //s //fast


//pump pulse:
#define laser_energy (6.04e3) //J/kg //Obergfell
#define laser_pulse_duration (50e-15) //s //Obergfell
//#define laser_pulse_duration (45e-15) //s
#define laser_wavelength (800e-9) //m //Obergfell
//#define laser_wavelength (800e-9) //m
#define laser_timepeak (0.0) //s

// rather irrelevant:
#define laser_fluence (34) //J/m² //Obergfell


//probe pulse:
#define probe_time_step_factor (500)
//#define probe_wavelength_amount_const (31)
#define probe_wavelength_amount_const (50)

#define min_probe_wavelength_const (1.25) //eV
//#define max_probe_wavelength_const (2.8) //eV
#define max_probe_wavelength_const (3.6) //eV
#define probe_wavelength_const (400e-9) //m //only for single probe wavelength use
#define laser_incidence_angle (7.0 * M_PI/180) //radian

#define probe_thickness_const (24e-9) //m //Obergfell
//#define probe_thickness_const (26.5e-9) //m //26.46e-9 // Z. Chen paper \pm10%

/*************************************************/

//2T constants:
#define interband_time_2T_const (10e-15) //s

//total absorption coefficient:
// Absorption coefficients vary too much, as to define them as a constant value
//#define absorption_coefficient (6.1355e7) //1/m // total absorption from https://refractiveindex.info/?shelf=main&book=Au&page=Johnson \hbar \omega = 3.1 eV
//#define absorption_coefficient (7.7089e7) //1/m // total absorption from https://refractiveindex.info/?shelf=main&book=Au&page=Johnson \hbar \omega = 1.55 eV

#endif // USED_DATA_H
