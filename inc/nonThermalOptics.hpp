

/*****************************
*File: nonThermalOptics.hpp  *
*Created on 01-03-2002       *
*Author: Ndione @ AG Rethfeld*			     	
******************************/



#ifndef NONTHERMALOPTICS_HHH
#define NONTHERMALOPTICS_HPP


#include <iostream>
#include <boost/filesystem.hpp>
#include <vector>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>



class nonThermalOptics{

	private:

		double *energy;
		double *distribution_f;
		
		gsl_interp_accel *f_accel;
    		gsl_spline  *f_spline;

		
	public:	
		nonThermalOptics();
		nonThermalOptics(std::string);
		~nonThermalOptics();
		
		std::vector<std::vector<double>> fileLoader;
		double interpolateDistribution(double);
	
};



#endif


