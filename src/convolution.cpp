
/************************************** 
 * File:   convolution.cpp
 * Author: Ndione, Ag Rethfeld
 **************************************
 * Created on July 02, 2019, 09:44 AM
 **************************************/

#include "convolution.h"

convolution::convolution(){

}
convolution::convolution(std::string filename){

	int_pointer = new Cintegration;

	if (filename != ""){

        	std::ifstream data(filename);

                std::vector <double>load(4);

                if (data.is_open()){

                        while(data >> load[0] >> load[1] >> load[2] >> load[3]){

                                if(LoadData.size()>0){

                                        load[0] *= 1.0;   //Time in fs and  

                                        load[1] *= 1.0;     //Reflectivity unitless

					load[2] *= 1.0;

					load[3] *= 1.0;
					
                                        LoadData.push_back(load);

                                }else{

                                        LoadData.push_back(load);

                                }

                        }
                }

        }else{
                std::cerr<<"\nThe file "<<filename<< " was not correctly opened or not found! Please check the file's name or it's path\n";
                std::exit(1);
        }


	this ->Time = new double[(int) LoadData.size()];

	this ->n_reflec  = new double[(int) LoadData.size()];

	this ->s_reflec = new double[(int) LoadData.size()];

	this ->p_reflec = new double[(int) LoadData.size()];
	
	for(int j = 0; j != (int) LoadData.size(); j++){

		Time[j] = LoadData[j][0] * 1e-15; // convert to second
			
		n_reflec[j] = LoadData[j][1];

		s_reflec[j] = LoadData[j][2];

		p_reflec[j] = LoadData[j][3];
	
		//std::cout<<Time[j]<<"\t"<<n_reflec[j]<<"\t"<<s_reflec[j]<<"\t"<<p_reflec[j]<<"\n";
		
	}
	


	n_reflec_accelerator  = gsl_interp_accel_alloc();

	n_reflec_spline = gsl_spline_alloc(gsl_interp_akima, (int) LoadData.size());

	gsl_spline_init(n_reflec_spline, Time, n_reflec, (int) LoadData.size());

	s_reflec_accelerator = gsl_interp_accel_alloc();

	s_reflec_spline = gsl_spline_alloc(gsl_interp_akima, (int) LoadData.size());

	gsl_spline_init(s_reflec_spline, Time, s_reflec, (int) LoadData.size());

	p_reflec_accelerator = gsl_interp_accel_alloc();

	p_reflec_spline = gsl_spline_alloc(gsl_interp_akima, (int) LoadData.size());

	gsl_spline_init(p_reflec_spline, Time, p_reflec, (int) LoadData.size());

	this ->lBound = Time[0];

	this ->uBound = Time[(int) LoadData.size() - 2];

	this ->FWHM = 540e-15;
}

convolution::~convolution(){

	delete int_pointer;

	delete [] Time;
	delete [] n_reflec;
	delete [] s_reflec;
	delete [] p_reflec;

	gsl_spline_free(n_reflec_spline);
	gsl_interp_accel_free(n_reflec_accelerator);
	gsl_spline_free(s_reflec_spline);
	gsl_interp_accel_free(s_reflec_accelerator);
	gsl_spline_free(p_reflec_spline);
	gsl_interp_accel_free(p_reflec_accelerator);		
	
}

double convolution::GetLowerBound(){

	return lBound;
}

double convolution::GetUpperBound(){

	return uBound;
}

double convolution::GetFWHM(){

	return FWHM;
}

double convolution::InterpReflectivity_n(double timing){

	if((int) LoadData.size() > 0){

		return gsl_spline_eval(n_reflec_spline, timing, n_reflec_accelerator);
	}else{

		std::cerr<<"\nThe file you wish to interpolate has no data points! The program exits\n";
		std::exit(1);	
	}	
}

double convolution::InterpReflectivity_sPol(double timing){

	if((int) LoadData.size() > 0){

		return gsl_spline_eval(s_reflec_spline, timing, s_reflec_accelerator);
	}else{

		std::cerr<<"\nThe file you wish to interpolate has no data points! The program exits\n";
		std::exit(1);	
	}	
}

double convolution::InterpReflectivity_pPol(double timing){

	if((int) LoadData.size() > 0){

		return gsl_spline_eval(p_reflec_spline, timing, p_reflec_accelerator);
	}else{

		std::cerr<<"\nThe file you wish to interpolate has no data points! The program exits\n";
		std::exit(1);	
	}	
}

double convolution::NormalizationFactor(){

	return std::sqrt((4.0*std::log(2))/(M_PI*GetFWHM()*GetFWHM()));
}

/********************************************************
The integration boundaries can be changed. One has
to care that the reflectivity file which is read 
covers the integration interval.
*********************************************************/
double convolution::ConvolutedReflectivity_n(double time, double time_initial, double time_final){
	 
	auto integrand = [ = ] (double x) ->double {

        	return  InterpReflectivity_n(x) *  exp(-4.0*log(2.0)*pow(((time-x)/GetFWHM()),2.0));
		 
	};

	auto integral = int_pointer->integrate_qag(integrand, time_initial , time_final);
	
	return NormalizationFactor() * integral;
}

double convolution::ConvolutedReflectivity_sPol(double time, double time_initial, double time_final){
	 
	auto integrand = [ = ] (double x) ->double {

        	return  InterpReflectivity_sPol(x) *  exp(-4.0*log(2.0)*pow(((time-x)/GetFWHM()),2.0));
		 
	};

	double integral = int_pointer->integrate_qag(integrand, time_initial , time_final);

	return NormalizationFactor() * integral;
}

double convolution::ConvolutedReflectivity_pPol(double time, double time_initial, double time_final){

	auto integrand = [ = ] (double x) ->double {

		return InterpReflectivity_pPol(x) * exp(-4.0*log(2.0)*pow(((time-x)/GetFWHM()),2.0));
	};

	double integral = int_pointer->integrate_qag(integrand, time_initial, time_final);
	 
	return NormalizationFactor() * integral;
}

double convolution::Derivative1_n(double time){

	if((int) LoadData.size() > 0){

		return gsl_spline_eval_deriv(n_reflec_spline, time, n_reflec_accelerator);
	}else{

		std::cerr<<"\nFirst derivative not found! The program exits\n";
		std::exit(1);	
	}		
}
		 
double convolution::Derivative2_n(double time){

	if((int) LoadData.size() > 0){

		return gsl_spline_eval_deriv2(n_reflec_spline, time, n_reflec_accelerator);
	}else{

		std::cerr<<"\nSecond derivative not found! The program exits\n";
		std::exit(1);	
	}	
}

double convolution::Derivative1_s(double time){

	if((int) LoadData.size() > 0){

		return gsl_spline_eval_deriv(s_reflec_spline, time, s_reflec_accelerator);
	}else{

		std::cerr<<"\nFirst derivative not found! The program exits\n";
		std::exit(1);	
	}
}

double convolution::Derivative2_s(double time){

	if((int) LoadData.size() > 0){

		return gsl_spline_eval_deriv2(s_reflec_spline, time, s_reflec_accelerator);
	}else{

		std::cerr<<"\nSecond derivative not found! The program exits\n";
		std::exit(1);	
	}
}

double convolution::Derivative1_p(double time){

	if((int) LoadData.size() > 0){

		return gsl_spline_eval_deriv(p_reflec_spline, time, p_reflec_accelerator);
	}else{

		std::cerr<<"\nFirst derivative not found! The program exits\n";
		std::exit(1);	
	}
}

double convolution::Derivative2_p(double time){

	if((int) LoadData.size() > 0){

		return gsl_spline_eval_deriv2(p_reflec_spline, time, p_reflec_accelerator);
	}else{

		std::cerr<<"\nSecond derivative not found! The program exits\n";
		std::exit(1);	
	}
}

double convolution::Integral_n(double lowerLimit, double upperLimit){

	if((int) LoadData.size() > 0){

		return gsl_spline_eval_integ(n_reflec_spline, lowerLimit, upperLimit, n_reflec_accelerator);
	}else{

		std::cerr<<"\nIntegral of interpolated function not found! The program exits\n";
		std::exit(1);	
	} 
}

double convolution::Integral_s(double lowerLimit, double upperLimit){

	if((int) LoadData.size() > 0){

		return gsl_spline_eval_integ(s_reflec_spline, lowerLimit, upperLimit, s_reflec_accelerator);
	}else{

		std::cerr<<"\nIntegral of interpolated function not found! The program exits\n";
		std::exit(1);	
	} 
}

double convolution::Integral_p(double lowerLimit, double upperLimit){

	if((int) LoadData.size() > 0){

		return gsl_spline_eval_integ(p_reflec_spline, lowerLimit, upperLimit, p_reflec_accelerator);
	}else{

		std::cerr<<"\nIntegral of interpolated function not found! The program exits\n";
		std::exit(1);	
	} 
}





















