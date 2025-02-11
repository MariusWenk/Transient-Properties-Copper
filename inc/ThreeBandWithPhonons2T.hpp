/* ************************************
 * File:   ThreeBandWithPhonons2T.hpp
 * Author: Ndione
 *
 * Created on June 17, 2020, 5:00 PM
 *************************************/

#ifndef THREEBANDWITHPHONONS2T_HPP
#define THREEBANDWITHPHONONS2T_HPP

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h> 
#include <boost/numeric/odeint.hpp>

#include "atom_unit_physics.h"
#include "dielectric_function.h"
#include "dielectric_function_drude.h"
#include "fermi_distribution.h"
#include "gsl_function_pp.h"
#include "integration.h"
#include "interpolation1d.hpp"
#include "interpolation2d.hpp"
#include "Laser.h"
#include "physical_constants.h"
#include "root_finding.h"
#include "ode.h"
#include "solver.h"


typedef std::vector<double> state_type;
typedef unsigned int uint;


class ThreeBandWithPhonons2T : public ODE{
    
	public:
		ThreeBandWithPhonons2T();    	
		ThreeBandWithPhonons2T(std::string sp_filename2D, std::string d_filename2D, Cfermi_distribution *fermi_ptr, uint col0, uint col1, uint col2); 
	    	virtual ~ThreeBandWithPhonons2T();

	    	
  
	private:
		
		
		
};




#endif /* THREEBANDWITHPHONONS2THPP */

