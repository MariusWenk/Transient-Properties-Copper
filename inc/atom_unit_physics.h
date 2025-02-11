/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   atom_unit_physics.h
 * Author: ndione, Ag Rethfeld
 *
 * Created on 5. Juni 2018, 16:21
 */

#ifndef ATOM_UNIT_PHYSICS_H
#define ATOM_UNIT_PHYSICS_H

#include "used_data.h"

//using namespace std;



#define fine_struct_const (7.2973525664e-3)
#define speed_of_light_au (1/fine_struct_const)
#define PI (3.14159265358979323846264338328)
#define m_e_au (1)
#define k_B_au (1)
//#define k_B_au (3.167e-6) /effectively this is implemented as K_in_au
#define hbar_au (1)
#define planck_const_au (hbar_au*2*PI)
#define vacuum_permittivity_au (1.0/(4.0*PI))


#define K_in_au (3.16681542254311e-6)
#define Joule_per_eV (1.602176487e-19)  
#define eV_to_Hartree (3.67493e-2)
#define Hartree_to_eV (27.2114) 
#define s_to_hbar_E_H (4.1341e16)
#define m_to_a0 (1.8898e10)
#define a0_to_m (5.2917e-11)
#define a3_to_m3 (a0_to_m*a0_to_m*a0_to_m)
#define n_au_to_SI (1.0/a3_to_m3)
#define Bohr_to_m (0.0529177211e-9)
#define Kg_in_au (1.0/9.10938188e-31)
#define speed_sound_au (3240*m_to_a0/s_to_hbar_E_H)
#define kg_in_au (3.2161e30)

#define effective_sp_mass_au (sp_band_mass_const*m_e_au)  
#define effective_d_mass_au (d_band_mass_const*m_e_au) 
#define unit_volume_au (unit_volume_SI*m_to_a0*m_to_a0*m_to_a0)

#define E_field_au_to_SI (5.14220674763e11)

#endif /* ATOM_UNIT_PHYSICS_H */

