

#ifndef  MAIN_H
#define MAIN_H


//#include <mpi.h>
#include <time.h>
#include<array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random> 
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/util/bind.hpp>
#include <boost/numeric/odeint/stepper/base/explicit_error_stepper_fsal_base.hpp>
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/default_operations.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/util/state_wrapper.hpp>
#include <boost/numeric/odeint/util/is_resizeable.hpp>
#include <boost/numeric/odeint/util/resizer.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54_classic.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>

#include "dielectric_function.h"
#include "dielectric_function_drude.h"
#include "dielectric_function_lindhard.h"
#include "dielectric_function_lorentz.h"
#include "dielectric_function_Mermin.h"
#include "fermi_distribution.h"
#include "transient_optical_param.h"
#include "sp_band.h"
#include "D_band.h"
#include "atom_unit_physics.h"
#include "physical_constants.h"
#include "ode.h"
#include "solver.h"
#include "transient_optical_param.h"
#include "BoostRootFinding1D.h"
#include "TwoBandWithPhonons.h"
#include "TwoBandWithPhonons2T.hpp"
#include "ThreeBandWithPhonons.h"
#include "ThreeBandWithPhonons2T.hpp"
#include "FourBandWithPhonons.h"
#include "FourBandWithPhonons2T.hpp"
#include "FourBandWithPhonons.h"
#include "distribution_neq.h"
#include "convolution.h"
#include "spline.h"
#include "nonThermalOptics.hpp"
#include "fitSilaeva.hpp"
#include "used_data.h"

 

 
#endif
