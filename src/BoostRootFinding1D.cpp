
 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/* 
 * File:   BoostRootFinding1D.cpp
 * Author: Ndione, Ag Rethfeld
 * 
 * Created on October 31, 2018, 1:52 PM
 */
#include "BoostRootFinding1D.h"

BoostRootFinding1D::BoostRootFinding1D(std::string dos_filename) {

	integrator_ptr = new Cintegration;

    	if (dos_filename != "") {

        	std::ifstream data(dos_filename);

        	std::vector<double> value(4);
        
		if (data.is_open()) {

            		while (data >> value[0] >> value[1] >> value[2] >> value[3]) {

                		if (value[3] != 0.) {

                    			if (dosfile.size() > 0) {
                        			
						value[0] -= dosfile[0][0]; //fermi energy is 0 in data
                        
						value[0] *= Hartree_to_eV;
                        
						value[1] /= Hartree_to_eV;
                        
						value[2] /= Hartree_to_eV;
                        
						value[3] /= Hartree_to_eV;
                        			dosfile.push_back(value);
                    			}else{
	
        			                dosfile.push_back(value);
                    			}
                		}
            		}
            		fermi_energy = (-1.0 * dosfile[0][0])  * Hartree_to_eV; // directly obtained from input Seb_DOS which had FermiEnergy at 0
            
			std::cout << "\nIn class \"BoostRootFinding1D\", Fermi energy located at " << fermi_energy << " eV in DOS file!\n";
            
			dosfile[0][0] = 0.0;
         
	   		data.close();
            		
			std::cout << "\nMax DOS energy located at " << dosfile[dosfile.size()-1][0] << " eV in total DOS file\n";
		
		}else{
            		std::cerr << "\nThe file: " << dos_filename << " for the realistic DOS is not correctly opened. \nPlease check the filename or it's path!\n";
            		std::exit(1);
        	}
         
	}
    //    for (unsigned int k = 0; k < (unsigned int) dosfile.size(); k++) {
    //
    //        std::cout<<"energy "<<dosfile[k][0]<<"\t"<<"dos_tot "<<dosfile[k][1]<<"\n";
    //         
    //    }
    
	this -> n = dosfile.size()-1;

    	this->energy = new double [(int) dosfile.size()];
    	
	this->DOS = new double [(int) dosfile.size()];
    	
	this->sp_DOS = new double [(int) dosfile.size()];
    	
	this->d_DOS = new double [(int) dosfile.size()];

	this ->distance_d_band_E_F = distance_d_band_E_F_const; // in eV

	this ->sp_band_mass = sp_band_mass_const * m_e;

	this ->d_band_mass = d_band_mass_const * m_e;	
 
	for (unsigned int i = 0; i < (unsigned int) dosfile.size(); i++) {
		energy[i] = dosfile[i][0];
		sp_DOS[i] = dosfile[i][1] * sp_DOS_stretching;
		d_DOS[i] = dosfile[i][2] * d_DOS_stretching; 
		//d_DOS[i] = dosfile[i][2] * SmoothFermiFunction(energy[i]) * d_DOS_stretching;  // d_DOS weighted  to obtain 10d-electrons  
		DOS[i] = sp_DOS[i] + d_DOS[i];
		
	}
	
	total_dos_accelerator = gsl_interp_accel_alloc();
	
	total_dos_spline = gsl_spline_alloc(gsl_interp_akima, (int) dosfile.size());
	gsl_spline_init(total_dos_spline, energy, DOS, (int) dosfile.size());
	
	sp_accel = gsl_interp_accel_alloc();
	sp_spline = gsl_spline_alloc(gsl_interp_akima, (int) dosfile.size());
	
	gsl_spline_init(sp_spline, energy, sp_DOS, (int) dosfile.size());
	d_accel = gsl_interp_accel_alloc();
	
	d_spline = gsl_spline_alloc(gsl_interp_akima, (int) dosfile.size());
	gsl_spline_init(d_spline, energy, d_DOS, (int) dosfile.size());

	this ->temperature_start = temperature_start_const;

	this ->LowerBound = energy[0];

	this ->UpperBound = energy[dosfile.size()-2]; 
   
	//this ->UpperBound = dosfile[dosfile.size() - 1] [0];

	this ->unitVolume = unit_volume_au / std::pow(m_to_a0, 3);
    
}

BoostRootFinding1D::BoostRootFinding1D(){
	this -> n = 4000;
	
	distance_d_band_E_F = distance_d_band_E_F_const;

	this ->sp_band_mass = sp_band_mass_const*m_e;

	this ->d_band_mass = d_band_mass_const*m_e;

	this ->fermi_energy = fermi_energy_const;

	this ->temperature_start = temperature_start_const;

	this ->unitVolume = unit_volume_au / std::pow(m_to_a0, 3);

	double delta_E = 2.0 * fermi_energy / n;

   	this->energy = new double [n];

        this->DOS = new double [n];

    	this->sp_DOS = new double [n];

    	this->d_DOS = new double [n];	
	 

	for (unsigned int i = 0; i < n; i++) {

        	energy[i] = i*delta_E;

        	sp_DOS[i] = sp_FEG(energy[i]*Joule_per_eV)*GetUnitVolume()*Joule_per_eV;

        	d_DOS[i] = d_FEG(energy[i]*Joule_per_eV)*GetUnitVolume()*Joule_per_eV;

        	DOS[i] = sp_DOS[i] + d_DOS[i];

		//std::cout<<energy[i]<<"\t"<<sp_DOS[i]<<"\t"<<d_DOS[i]<<"\t"<<DOS[i]<<"\n";  
	} 
}

BoostRootFinding1D::~BoostRootFinding1D() {
	 
    	delete [] energy;
    	delete [] sp_DOS;
    	delete [] d_DOS;
    	delete [] DOS;
    	delete integrator_ptr;
	 
	gsl_spline_free(total_dos_spline);
	gsl_interp_accel_free(total_dos_accelerator);
	gsl_spline_free(sp_spline);
	gsl_interp_accel_free(sp_accel);
	gsl_spline_free(d_spline);
	gsl_interp_accel_free(d_accel);
}

double BoostRootFinding1D::get_sp_bandMass(){

	return sp_band_mass;
}

double BoostRootFinding1D::get_d_bandMass(){

	return d_band_mass;
}

double BoostRootFinding1D::GetUnitVolume(){

	return unitVolume;
}

double BoostRootFinding1D::GetDistance_d_BandtoFermi(){

	return distance_d_band_E_F;
}

double BoostRootFinding1D::GetFermiEnergy() {

    	return fermi_energy;
}

double BoostRootFinding1D::GetTemperatureStart(){

	return temperature_start;
}

double BoostRootFinding1D::GetLowerBound() {

    	return LowerBound;
}

double BoostRootFinding1D::GetUpperBound() {

    	return UpperBound;
}

double BoostRootFinding1D::GetNumberDiscPoints(){

	return n;
}

double BoostRootFinding1D::GetFermiDistribution(double energy, double T, double mu) {
    
	return 1.0 / (std::exp((energy - mu) / (k_B_eV * T)) + 1.0);
	//	return 1.0 / (std::exp((energy - mu) / (k_B * T)) + 1.0);
}

double BoostRootFinding1D::StepFunction(double energy){
	
	if ((energy >= 0.0) && (energy < (GetFermiEnergy() -  GetDistance_d_BandtoFermi()))){

		return  1.0;
	} 	
	return 0.0;	
}

double BoostRootFinding1D::SmoothFermiFunction(double energy){
	
	double temp = temperature_smooth_fermi_function;
		
	//double mu = GetFermiEnergy() -  GetDistance_d_BandtoFermi();
	double mu = GetFermiEnergy() -  1;	 

	return GetFermiDistribution(energy, temp, mu);
}

double BoostRootFinding1D::sp_FEG(double energy){

	return (std::sqrt(2.0 * (energy - 0.0)) * std::pow(get_sp_bandMass(), 1.5)) / (PI * PI * std::pow(hbar, 3.0)); 	
}

double BoostRootFinding1D::d_FEG(double energy){
 
	if ((energy >= 0) && (energy < ((GetFermiEnergy()*Joule_per_eV)-(GetDistance_d_BandtoFermi()*Joule_per_eV)))) {
	
        	return (std::sqrt(2.0 * (energy - 0.0)) * std::pow(get_d_bandMass(), 1.5)) / (PI * PI * std::pow(hbar, 3.0));
    	}
         
	return 0.0;	  
}

double BoostRootFinding1D::GetDOS(double energy) {

	if (dosfile.size() > 0) {
        	
		if ((energy < dosfile[0][0]) || (energy > dosfile[dosfile.size() - 2][0])) {
   
           		 std::cerr<<"\nEnergy range not supported by the interpolation! The program exits!\n";
				 std::exit(1);
        	}else{
	
			return gsl_spline_eval(total_dos_spline, energy, total_dos_accelerator);
		}
        
	}else{
       		 std::cerr << "\nNo DOS specified or DOS doesn't contain any discret point!\n"
                 << "The program exits!\n";
       		 std::exit(1);
	}
}

double BoostRootFinding1D::GetDOS_sp(double energy) {

	if (dosfile.size() > 0) {
        	
		if ((energy < dosfile[0][0]) || (energy > dosfile[dosfile.size() - 2][0])) {
   
           		 std::cerr<<"\nEnergy range not supported by the interpolation! The program exits!\n";
				 std::exit(1);
        	}else{
	
			return gsl_spline_eval(sp_spline, energy, sp_accel);
		}
        
	}else{
       		 std::cerr << "\nNo sp_DOS specified or DOS doesn't contain any discret point!\n"
                 << "The program exits!\n";
       		 std::exit(1);
	}
}
double BoostRootFinding1D::GetDOS_d(double energy) {

	if (dosfile.size() > 0) {
        	
		if ((energy < dosfile[0][0]) || (energy > dosfile[dosfile.size() - 2][0])) {
   
           		 std::cerr<<"\nEnergy range not supported by the interpolation! The program exits!\n";
				 std::exit(1);
        	}else{
	
			return gsl_spline_eval(d_spline, energy, d_accel) ;
		}
        
	}else{
       		 std::cerr << "\nNo d_DOS specified or DOS doesn't contain any discret point!\n"
                 << "The program exits!\n";
       		 std::exit(1);
	}
}
double BoostRootFinding1D::NumberofParticles(double T, double mu) {
	
	/*
	auto integrand = [ = ](double energy) -> double {
        
		return GetDOS(energy)*GetFermiDistribution(energy, T, mu);				
         };
        
	        //return integrator_ptr->integrate_qag(integrand, GetLowerBound(), GetUpperBound());
		return integrator_ptr->integrate_qag(integrand, 0.0, GetFermiEnergy());
					 
    	*/	
   	double sum = 0.0;

    	if (T > 0.0) {

        	for (unsigned int i = 1; i < n - 1; i++) {

			sum += (DOS[i] * GetFermiDistribution(energy[i], T, mu))*(energy[i + 1] - energy[i - 1]);
        	}
        
		return 0.5 * (sum
	                + DOS[0] * GetFermiDistribution(energy[0], T, mu)*(energy[1] - energy[0])
	                + DOS[n - 1] * GetFermiDistribution(energy[n - 1], T, mu)*(energy[n - 1] - energy[n - 2])
	                );
    	}else{

	        for (unsigned int i = 1; i < n - 1; i++) {
	
	            	sum += (energy[i] > fermi_energy) ? 0.0 : DOS[i]*(energy[i + 1] - energy[i - 1]);
        	}

        	return 0.5 * (sum + DOS[0]*(energy[1] - energy[0]) + DOS[n - 1]*(energy[n - 1] - energy[n - 2]));
    	}
}
double BoostRootFinding1D::sp_NumberofParticles(double T, double mu) {
	
	/*auto integrand = [ = ](double energy) -> double {
        
		return GetDOS_sp(energy)*GetFermiDistribution(energy, T, mu);
	};
        
	//return integrator_ptr->integrate_qag(integrand, GetLowerBound(), GetUpperBound());
	return integrator_ptr->integrate_qag(integrand, 0, GetFermiEnergy());
   	
	*/
   	double sum = 0.0;

    	if (T > 0.0) {

	        for (unsigned int i = 1; i < n - 1; i++) {
            	
			sum += (sp_DOS[i] * GetFermiDistribution(energy[i], T, mu))*(energy[i + 1] - energy[i - 1]);
        	}
        
		return 0.5 * (sum
                	+ sp_DOS[0] * GetFermiDistribution(energy[0], T, mu)*(energy[1] - energy[0])
                	+ sp_DOS[n - 1] * GetFermiDistribution(energy[n - 1], T, mu)*(energy[n - 1] - energy[n - 2])
                	);
    	}else{

        	for (unsigned int i = 1; i < n - 1; i++) {

            		sum += (energy[i] > fermi_energy) ? 0.0 : sp_DOS[i]*(energy[i + 1] - energy[i - 1]);
        	}
        
		return 0.5 * (sum + sp_DOS[0]*(energy[1] - energy[0]) + sp_DOS[n - 1]*(energy[n - 1] - energy[n - 2]));
    	}
}

double BoostRootFinding1D::d_NumberofParticles(double T, double mu) {

	/*
           auto integrand = [ = ](double energy) -> double {
        
                    return GetDOS_d(energy)*GetFermiDistribution(energy, T, mu);
                };
        
	        //return integrator_ptr->integrate_qag(integrand, GetLowerBound(), GetUpperBound());
		return integrator_ptr->integrate_qag(integrand, 0, GetFermiEnergy());
      	
	*/
   	double sum = 0.0;

    	if (T > 0.0) {

        	for (unsigned int i = 1; i < n - 1; i++) {

            		sum += (d_DOS[i] * GetFermiDistribution(energy[i], T, mu))*(energy[i + 1] - energy[i - 1]);
        	}
        
		return 0.5 * (sum
                	+ d_DOS[0] * GetFermiDistribution(energy[0], T, mu)*(energy[1] - energy[0])
                	+ d_DOS[n - 1] * GetFermiDistribution(energy[n - 1], T, mu)*(energy[n - 1] - energy[n - 2])
                	);
    
	}else{

        	for (unsigned int i = 1; i < n - 1; i++) {
            		
			sum += (energy[i] > fermi_energy) ? 0.0 : d_DOS[i]*(energy[i + 1] - energy[i - 1]);
        	}
        
		return 0.5 * (sum + d_DOS[0]*(energy[1] - energy[0]) + d_DOS[n - 1]*(energy[n - 1] - energy[n - 2]));
    	}
}

double BoostRootFinding1D::df_dmu(double energy, double T, double mu) {

    	return (1.0 / (k_B_eV * T)) / (2.0 + 2.0 * std::cosh((energy - mu) / (k_B_eV * T)));
	//return (1.0 / (k_B * T)) / (2.0 + 2.0 * std::cosh((energy - mu) / (k_B * T)));
}

double BoostRootFinding1D::Get_P_mu(double T, double mu) {

   	/*  if (T > 0.) {
        auto integrand = [ = ] (double energy) ->double {
            return GetDOS(energy) * df_dmu(energy, T, mu);
        };
        return integrator_ptr->integrate_qag(integrand, GetLowerBound(), GetUpperBound());
    	} else {
        return 0.0;
    	}*/
	double sum = 0.;

     	if(T>0.0){	       		
    
        	for (unsigned int i = 1; i < n - 1; i++) {
    
			sum += (DOS[i] / (2.0 + 2.0 * std::cosh((energy[i] - mu) / (k_B_eV * T)))*(energy[i + 1] - energy[i - 1]));
        	}
    
        	return 0.5 / (k_B_eV * T)*
                	(
                	sum
                	+ DOS[0] / (2.0 + 2.0 * std::cosh((energy[0] - mu) / (k_B_eV * T)))*(energy[1] - energy[0])
                	+ DOS[n - 1] / (2.0 + 2 * std::cosh((energy[n - 1] - mu) / (k_B_eV * T)))*(energy[n - 1] - energy[n - 2])
                	);
     	}else{
	
		return 0.0;
    	} 
		
	/*if(T>0.0){	
       double sum = 0.;
    
        for (unsigned int i = 1; i < n - 1; i++) {
    
    
            sum += (DOS[i] / (2.0 + 2.0 * std::cosh((energy[i] - mu) / (k_B * T)))*(energy[i + 1] - energy[i - 1]));
        }
    
        return 0.5 / (k_B * T)*
                (
                sum
                + DOS[0] / (2.0 + 2.0 * std::cosh((energy[0] - mu) / (k_B * T)))*(energy[1] - energy[0])
                + DOS[n - 1] / (2.0 + 2 * std::cosh((energy[n - 1] - mu) / (k_B * T)))*(energy[n - 1] - energy[n - 2])
                );
     }else {
	return 0.0;
    }*/
}

/*******************************************
*The chemical potential is provided in ev
*It needs as argument temperature in Kelvin
********************************************/
double BoostRootFinding1D::GetMuRootFinding(double T) {

    	double N0 = NumberofParticles(GetTemperatureStart(), GetFermiEnergy());
	
	using namespace boost::math::tools;

        double guess = GetFermiEnergy(); //initial guess

        double min = 0.75 * GetFermiEnergy(); //Minimum possible value

        double max = 10.0 * GetFermiEnergy(); //Maximum possible value

        const int digits = std::numeric_limits<double>::digits; //Maximum possible binary digits accuracy for type "double"

        int get_digits = static_cast<int> (digits * 0.6); //Accuracy doubles with each step, so stop when ever half of the digits are correct

        const boost::uintmax_t maxit = 20; //Maximum number of iterations

        boost::uintmax_t it = maxit; //Initally our chosen max iterations, but updated with actual.

	if (T > 0.) {

	        double result = newton_raphson_iterate([ = ](const double& mu){

            	return std::make_pair(NumberofParticles(T, mu) - N0, Get_P_mu(T, mu));

        	}, guess, min, max, get_digits, it);

	        return result;
    	}

        return GetFermiEnergy();
    	
}
