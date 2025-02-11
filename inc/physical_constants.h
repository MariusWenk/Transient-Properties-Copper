#ifndef PHYSICAL_CONSTANTS_H
#define PHYSICAL_CONSTANTS_H

#include "used_data.h"


// Integration parameters 
//#define interval (100000)
//#define estimated_integration_error (-1)
//#define integration_routine (GSL_INTEG_GAUSS61)
//#define absolute_error (0)
//#define relative_error (1.0e-9)
//#define workspace_allocation (gsl_integration_workspace_alloc (interval))

#define speed_of_light (299792458) // m s^-1
#define PI (3.14159265358979323846264338328)
#define avogadro_const (6.02214199e23)
#define m_e (9.10938188e-31) // kg
#define electron_mass (m_e) // Kg
#define electron_effective_mass (1.0*electron_mass) //Kg // Can be modified according to the material which is used
#define electron_charge (-1.60217733e-19) // C
#define boltzmann_const (1.38064852e-23)  // J K^-1
#define k_B (boltzmann_const)
#define k_B_eV (8.617330e-5) //eV/K
#define vacuum_permittivity (8.854187817e-12) // F m^-1
#define planck_constant (6.62607004e-34)  // m^2 kg s^-1
#define reduced_planck_constant (planck_constant/(2.0*PI))  // m^2 kg s^-1
//#define temperature (3000.0) // K
#define vacuum_impedance (377) // Ohm
#define Joule_per_eV (1.602176487e-19) /* kg m^2 / s^2 */
//#define Hartree_to_eV (27.22138602)  // 1E_h=27.....eV
#define Joule_to_Hartree ( 1.0/4.35974417e-18) //1J=2.....e17
#define Hartree_to_Joule (4.3597443419e-18) // 1E_h = 3e.....e-18J
#define Kelvin_per_eV (11604.5250061657) // Kelivin K
#define hbar (1.054571800e-34) /* kg m^2 / s */
#define hbar_eV (6.5821195146e-16)
#define Hz_to_eV (4.13566553853599e-15)

#define SI_hbar (1.05457162825e-34) /* kg m^2 / s */
#define SI_k_B (1.3806504e-23) /* kg m^2 / K s^2 */
#define SI_speed_of_light (2.99792458e8) /* m / s */
#define SI_Joule_per_eV (1.602176487e-19) /* kg m^2 / s^2 */
#define SI_avogadro (6.02214199e23) /* 1 / mol */
#define SI_charge (1.602176487e-19) /* A s */
#define SI_vac_perm (8.854187817e-12)  /* A^2 s^4 / kg m^3 */
#define SI_m_e (9.10938188e-31) /*kg*/
#define SI_Hartree (4.359744650e-18)
#define SI_Bohr (0.52917721067e-10)
#define SI_temperature (31577513e-2) /* 1 Hartree in K*/
#define SI_eV_per_Hartree (SI_Hartree/SI_Joule_per_eV) //ca 27.2eV

#define AU_hbar (1) //hbar
#define AU_k_B (1) //kB
#define AU_speed_of_light (1/7.2973525664e-3) //E_H*a_0 / hbar
#define AU_m_e (1) //m_e
#define AU_charge (1) //e
#define AU_time (SI_Hartree/SI_hbar) // hbar/E_H, 1 second in AU
#define AU_density (1.481847107e-31) // = a_0^3, 1 m^-3 in AU

#define LOG_OF_INF (700.0)
#define std_temp (temperature_start_const/SI_temperature)

#endif /* PHYSICAL_CONSTANTS_H */

