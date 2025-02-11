
/************************************** 
 * File:   BoostRootFinding1D.h
 * Author: Ndione, Ag Rethfeld
 **************************************
 * Created on July 02, 2019, 09:44 AM
 **************************************/



#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <functional>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>


#include "integration.h"
#include "physical_constants.h"
#include "atom_unit_physics.h"

class convolution{

	public:
		convolution();
		convolution(std::string filename);
		virtual  ~convolution();

		//vector of vector containing information on time and reflectivity
		std::vector<std::vector<double>> LoadData;
		//! Hands over integration bounds
		double GetLowerBound();
		double GetUpperBound();
		//Interpolation function for normal incidence reflectivity
		double InterpReflectivity_n(double);
		//First and second derivative of an interpolated function
		double Derivative1_n(double);
		double Derivative2_n(double); 
		//Interpolation function for s-polarized reflectivity
		double InterpReflectivity_sPol(double);
		//First and second derivative of an interpolated function
		double Derivative1_s(double);
		double Derivative2_s(double);
		//Interpolation function for p-polarized reflectivity
		double InterpReflectivity_pPol(double);
		//First and second derivative of an interpolated function
		double Derivative1_p(double);
		double Derivative2_p(double);
		//!Hands over the numerical integral of an interpolated function over a range [a, b]
		double Integral_n(double, double);
		double Integral_s(double, double);
		double Integral_p(double, double);

		//!Hands over the convolution of reflectivity signal a Gaussian
		double ConvolutedReflectivity_n(double, double, double);
		double ConvolutedReflectivity_sPol(double, double, double);
		double ConvolutedReflectivity_pPol(double, double, double);
		double GetFWHM();
		double NormalizationFactor();
 
	private:

		//gsl accelarators
		gsl_interp_accel *n_reflec_accelerator;
		gsl_interp_accel *s_reflec_accelerator;	
		gsl_interp_accel *p_reflec_accelerator;
		//gsl splines
		gsl_spline *n_reflec_spline;
		gsl_spline *s_reflec_spline;
		gsl_spline *p_reflec_spline;

		//Arrays for time and transient reflectivity
		double *Time;
		double *n_reflec;
		double *s_reflec;
		double *p_reflec;
		//Integration pointer
		Cintegration  *int_pointer;
		//Integration Bounds
		double lBound;
		double uBound;

		//FWHM for Gaussian
		double FWHM;
};		



#endif
