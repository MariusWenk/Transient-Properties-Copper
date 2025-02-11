

 
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h> 


#include "integration.h"
#include "interpolation1d.hpp"
#include "dielectric_function.h"
#include "physical_constants.h"
#include "gsl_function_pp.h"
#include "integration.h"
#include "atom_unit_physics.h"
#include "BoostRootFinding1D.h"
#ifndef SP_BAND_H
#define SP_BAND_H

class Sp_band {

public:

    //! Constructor
    Sp_band();
    //! Copy constructor
    Sp_band(int number_disc_points);
    //    Sp_band(const sp_band& orig);
    //! Destructor
    virtual ~Sp_band();
    
    //! Hands over s density after reading a total DOS file
    void Readfile_s_dos(std::vector<double> *energy, std::vector<double> *s_dos, const char *filename);
    //Hands over p density after reading a total DOS file
    void Readfile_p_dos(std::vector<double> *energy, std::vector<double> *p_dos, const char *filename);
    //Hands over sp density after reading a total DOS file
    void Readfile_sp_dos(std::vector<double> *energy, std::vector<double> *sp_dos, const char *filename);
    
    //! Hands over the FEG DOS at a given energy
    double Get_sp_DOS(double energy);
    //! Hands over the FEG DOS at a given number of discretization points
    double Get_sp_DOS(int number_of_disc_points);
    //! Hands over the FEG DOS at a given wavevector 
    double Get_sp_DOSWavevector(double wavevector);
    //! Hands over the band mass
    double Get_sp_BandMass();
    //! Hands over the unit cell volume of the material
    double GetUnitVolume();
    //! Hands over electron density of sp band
    double Get_n_sp_Density();
    //! hands over the number of discretization points
    int GetNumberDiscPoints();
    //! Hands over fermi energy of the material
    double GetFermiEnergy();
    //! Hands over the dispersion relation of a FEG
    double GetDispersionRelation(double wavevector);
    //! Hands over the wavevector
    double GetWavevector(double energy);
    void Read();

private:

    // Effective mass of sp band
    double sp_band_mass;
    // Number of discretization points
    int number_disc_points;
    // Unit cell volume
    double unit_cell_volume;
    // Free electron density
    double n_sp_density;
    // Fermi energy
    double fermi_energy;

};

#endif /* SP_BAND_H */

