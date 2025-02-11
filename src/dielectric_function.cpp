/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
   
/**************************************** 
 * File:   dielectric_function.cpp
 * Author: ndione, Ag Rethfeld
 * 
 * Created on February 9, 2017, 11:17 AM
 ***************************************/

#include "dielectric_function.h"


dielectric_function::dielectric_function(){
	
	integrator = new Cintegration;

	this ->angle_incidence = laser_incidence_angle * M_PI / 180.0; // convert degree to radian

	this ->film_thickness = probe_thickness_const;  //Usually measurement instruments for film thickness have 10% precision
} 

dielectric_function::~dielectric_function() {
    
    delete integrator;
}

// Functions have to be defined but should not be used(Should not be called from this class)

double dielectric_function::GetRealPart(double omega, double Te, double Tph) {
   
    	std::cerr <<"\nTried to call non-existing function\n"
            << "dielectric_function::GetRealPart(double omega, double Te, double Tph)\n"
            << "with omega = " << omega << "\n"
	    << "with Te = " << Te << "\n"
	    << "with Tph = " << Tph << "\n"
            << "Please use one of the derived classes.\n"
            << "Exiting program.\n";
    	exit(1);
}

double dielectric_function::GetRealPart(double omega, double Te) {
    
    	std::cerr << "\nTried to call non-existing function GetRealPart(double omega, double Te)\n"
            << "dielectric_function::\n"
            << "with omega = " << omega << "\n"
	    << "with Te = " << Te << "\n"
            << "Please use one of the derived classes.\n"
            << "Exiting program.\n";
    	exit(1);
}

double dielectric_function::GetRealPart(double omega, double Te, double Tph, double wavevector) {
   
	std::cerr << "\nTried to call non-existing function \n"
            << "dielectric_function::GetRealPart(double omega, double Te, double Tph, double wavevector)\n"
            << "with omega = " << omega << "\n"
	    << "with Te = " << Te << "\n"
	    << "with Tph = " << Tph << "\n"
            << "with wavevector = " << wavevector << "\n"
            << "Please use one of the derived classes.\n"
            << "Exiting program.\n";
    	exit(1);
}

double dielectric_function::GetImaginaryPart(double omega, double Te, double Tph) {
    
    	std::cerr <<"\nTried to call non-existing function\n"
            << "dielectric_function::GetImaginaryPart(double omega, double Te, double Tph)\n"
            << "with omega = " << omega << "\n"
	    << "with Te = " << Te << "\n"
	    << "with Tph = " << Tph << "\n"
            << "Please use one of the derived classes.\n"
            << "Exiting program.\n";
    	exit(1);
}

double dielectric_function::GetImaginaryPart(double omega, double Te){
 
    	std::cerr << "\nTried to call non-existing function \n"
            << "dielectric_function::GetImaginaryPart(double omega, double Te)\n"
            << "with omega = " << omega << "\n"
	    << "with Te = " << Te << "\n"
            << "Please use one of the derived classes.\n"
            << "Exiting program.\n";
    	exit(1);
}
 
double dielectric_function::GetImaginaryPart(double omega, double Te, double Tph, double wavevector){
 
	std::cerr << "\nTried to call non-existing function \n"
            << "dielectric_function::GetImaginaryPart(double omega, double Te, double Tph, double wavevector)\n"
            << "with omega = " << omega << "\n"
	    << "with Te = " << Te << "\n"
	    << "with Tph = " << Tph << "\n"
            << "with wavevector = " << wavevector << "\n"
            << "Please use one of the derived classes.\n"
            << "Exiting program.\n";
    	exit(1);
}

double dielectric_function::GetRealFromImaginaryPart(double omega, double Te, double Tph){
 
 	std::cerr <<"\nTried to call non-existing function\n"
            << "dielectric_function::GetRealFromImaginaryPart(double omega, double Te, double Tph)\n"
            << "with omega = " << omega << "\n"
	    << "with Te = " << Te << "\n"
	    << "with Tph = " << Tph << "\n"
            << "Please use one of the derived classes.\n"
            << "Exiting program.\n";
    	exit(1);
}

double dielectric_function::GetRealFromImaginaryPart(double omega, double Te){
 
    	std::cerr << "\nTried to call non-existing function \n"
            << "dielectric_function::GetRealFromImaginaryPart(double omega, double Te)\n"
            << "with omega = " << omega << "\n"
	    << "with Te = " << Te << "\n"
            << "Please use one of the derived classes.\n"
            << "Exiting program.\n";
    	exit(1);
}

double dielectric_function::GetRealFromImaginaryPart(double omega, double Te, double Tph, double wavevector){
 
    	std::cerr << "\nTried to call non-existing function \n"
            << "dielectric_function::GetRealFromImaginaryPart(double omega, double Te, double Tph, double wavevector)\n"
            << "with omega = " << omega << "\n"
	    << "with Te = " << Te << "\n"
	    << "with Tph = " << Tph << "\n"
            << "with wavevector = " << wavevector << "\n"
            << "Please use one of the derived classes.\n"
            << "Exiting program.\n";
    	exit(1);
}

double dielectric_function::GetAngleIncidence(){
	
	return angle_incidence;
}

double dielectric_function::GetFilmThickness(){

	return film_thickness;
}

std::complex<double> dielectric_function::n_top(std::complex<double> epsilon){

	return std::complex<double>(1.0,0.0); // Assuming vacuum. if not vacuum should return sqrt(epsilon)
}

std::complex<double> dielectric_function::n_bottom(std::complex<double> epsilon){

	return std::complex<double>(1.0,0.0); // Assuming vacuum. if not vacuum should return sqrt(epsilon)
}

std::complex<double> dielectric_function::n_film(std::complex<double> epsilon){

	return std::sqrt(epsilon);
}

double dielectric_function::TransmittedAngle(){

	return GetAngleIncidence();
}

std::complex<double>dielectric_function::cosineRefractedAngle(std::complex<double> epsilon){

	return sqrt(1.0-(pow((sin(GetAngleIncidence())/n_film(epsilon)),2.0)));
}

std::complex<double> dielectric_function::BetaFunction(double omega, std::complex<double> epsilon){

	//convert omega into wavelength
	double probe_wavelength = 2.0*M_PI*hbar_eV*speed_of_light/omega;  
	return (2.0*M_PI*GetFilmThickness()*n_film(epsilon)*cosineRefractedAngle(epsilon)) / probe_wavelength;

}

std::complex<double> dielectric_function::r_top_film_normalIncidence(std::complex<double> epsilon){

	std::complex<double> num = n_top(epsilon) - n_film(epsilon) ;
	std::complex<double> denom = n_top(epsilon) + n_film(epsilon) ;
	return num / denom;
}
std::complex<double> dielectric_function::r_top_film_sPolarized(std::complex<double> epsilon){

	std::complex<double> num = n_top(epsilon)*cos(GetAngleIncidence()) - n_film(epsilon)*cosineRefractedAngle(epsilon);
	std::complex<double> denom = n_top(epsilon)*cos(GetAngleIncidence()) + n_film(epsilon)*cosineRefractedAngle(epsilon);
	return num / denom;
}

std::complex<double> dielectric_function::t_top_film_sPolarized(std::complex<double> epsilon){

	std::complex<double> num = 2.0*cos(GetAngleIncidence());
	std::complex<double> denom =cos(GetAngleIncidence()) + n_film(epsilon)*cosineRefractedAngle(epsilon) ;
	return num / denom;
}

std::complex<double> dielectric_function::r_top_film_pPolarized(std::complex<double> epsilon){

	std::complex<double> num = n_film(epsilon)*cos(GetAngleIncidence()) - n_top(epsilon)*cosineRefractedAngle(epsilon);
	std::complex<double> denom = n_film(epsilon)*cos(GetAngleIncidence()) + n_top(epsilon)*cosineRefractedAngle(epsilon);
	return num / denom;
}

std::complex<double> dielectric_function::t_top_film_pPolarized(std::complex<double> epsilon){

	std::complex<double> num = 2.0*cos(GetAngleIncidence());
	std::complex<double> denom = n_film(epsilon)*cos(GetAngleIncidence()) + cosineRefractedAngle(epsilon);
	return num / denom;
}

std::complex<double> dielectric_function::t_top_film_normalIncidence(std::complex<double> epsilon){

	std::complex<double> num = 2.0*n_top(epsilon); // n_top is 1 for vacuum
	std::complex<double> denom = n_film(epsilon)  + n_top(epsilon);
	return num / denom;
}

std::complex<double> dielectric_function::r_film_bottom_normalIncidence(std::complex<double> epsilon){

	std::complex<double> num = n_film(epsilon) - n_bottom(epsilon);
	std::complex<double> denom = n_film(epsilon) + n_bottom(epsilon);
	return num / denom;
}

std::complex<double> dielectric_function::r_film_bottom_sPolarized(std::complex<double> epsilon){

	std::complex<double> num = n_film(epsilon)*cosineRefractedAngle(epsilon) - n_bottom(epsilon)*cos(TransmittedAngle());
	std::complex<double> denom = n_film(epsilon)*cosineRefractedAngle(epsilon) + n_bottom(epsilon)*cos(TransmittedAngle());
	return num / denom;
}

std::complex<double> dielectric_function::t_film_bottom_sPolarized(std::complex<double> epsilon){

	std::complex<double> num = 2.0*n_film(epsilon)*cosineRefractedAngle(epsilon);
	std::complex<double> denom = n_film(epsilon)*cosineRefractedAngle(epsilon) + cos(TransmittedAngle());
	return num / denom;
}

std::complex<double> dielectric_function::r_film_bottom_pPolarized(std::complex<double> epsilon){

	std::complex<double> num = n_bottom(epsilon)*cosineRefractedAngle(epsilon) - n_film(epsilon)*cos(TransmittedAngle());
	std::complex<double> denom = n_bottom(epsilon)*cosineRefractedAngle(epsilon) + n_film(epsilon)*cos(TransmittedAngle());
	return num / denom;
}

std::complex<double> dielectric_function::t_film_bottom_pPolarized(std::complex<double> epsilon){

	std::complex<double> num = 2.0*n_film(epsilon)*cosineRefractedAngle(epsilon);
	std::complex<double> denom = cosineRefractedAngle(epsilon) + n_film(epsilon)*cos(TransmittedAngle());
	return num / denom;
}

std::complex<double> dielectric_function::t_film_bottom_normalIncidence(std::complex<double> epsilon){

	std::complex<double> num = 2.0*n_film(epsilon)*cosineRefractedAngle(epsilon);
	std::complex<double> denom = n_bottom(epsilon) + n_film(epsilon);
	return num / denom;
}

double dielectric_function::sPolarizedReflectivityFilm(double omega, std::complex<double> epsilon){

	std::complex<double> complex_exp (cos(2.0*BetaFunction(omega, epsilon).real())* exp(-2.0*BetaFunction(omega, epsilon).imag()), sin(2.0*BetaFunction(omega, epsilon).real())* exp(-2.0*BetaFunction(omega, epsilon).imag()));
 
 
	return std::norm(		
		(r_top_film_sPolarized(epsilon)+ r_film_bottom_sPolarized(epsilon)*complex_exp)
		/
		(1.0 + r_top_film_sPolarized(epsilon)*r_film_bottom_sPolarized(epsilon)*complex_exp)
	);;
}

double dielectric_function::pPolarizedReflectivityFilm(double omega, std::complex<double> epsilon){

	std::complex<double> complex_exp (cos(2.0*BetaFunction(omega, epsilon).real())* exp(-2.0*BetaFunction(omega, epsilon).imag()), sin(2.0*BetaFunction(omega, epsilon).real())* exp(-2.0*BetaFunction(omega, epsilon).imag()));
 
 
	return std::norm(		
		(r_top_film_pPolarized(epsilon)+ r_film_bottom_pPolarized(epsilon)*complex_exp)
		/
		(1.0 + r_top_film_pPolarized(epsilon)*r_film_bottom_pPolarized(epsilon)*complex_exp)
	);
 
}

double dielectric_function::NormalIncidenceReflectivityFilm(double omega, std::complex<double> epsilon){

		std::complex<double> complex_exp (cos(2.0*BetaFunction(omega, epsilon).real())* exp(-2.0*BetaFunction(omega, epsilon).imag()), sin(2.0*BetaFunction(omega, epsilon).real())* exp(-2.0*BetaFunction(omega, epsilon).imag()));
 
 
	return std::norm(		
		(r_top_film_normalIncidence(epsilon)+ r_film_bottom_normalIncidence(epsilon)*complex_exp)
		/
		(1.0 + r_top_film_normalIncidence(epsilon)*r_film_bottom_normalIncidence(epsilon)*complex_exp)
	);
}

double dielectric_function::NormalIncidenceReflectivityFilm(double omega, double Te, double Tph, std::complex<double> epsilon){

	return NormalIncidenceReflectivityFilm(omega, epsilon);
}

double dielectric_function::pPolarizedReflectivityFilm(double omega, double Te, double Tph, std::complex<double> epsilon){

	return pPolarizedReflectivityFilm(omega, epsilon);
}

double dielectric_function::sPolarizedReflectivityFilm(double omega, double Te, double Tph, std::complex<double> epsilon){

	return sPolarizedReflectivityFilm(omega, epsilon);
}

double dielectric_function::sPolarizedTransmissivityFilm(double omega, std::complex<double> epsilon){

	std::complex<double> complex_exp2 (cos(2.0*BetaFunction(omega, epsilon).real())* exp(-2.0*BetaFunction(omega, epsilon).imag()), sin(2.0*BetaFunction(omega, epsilon).real())* exp(-2.0*BetaFunction(omega, epsilon).imag()));

	std::complex<double> complex_exp1 (cos(BetaFunction(omega, epsilon).real())* exp(-BetaFunction(omega, epsilon).imag()), sin(BetaFunction(omega, epsilon).real())* exp(-BetaFunction(omega, epsilon).imag()));

	return std::norm(
			(t_top_film_sPolarized(epsilon)*t_film_bottom_sPolarized(epsilon)*complex_exp1)
			/
			(1.0 + r_top_film_sPolarized(epsilon)*r_film_bottom_sPolarized(epsilon)*complex_exp2)
		);
	
}

double dielectric_function::pPolarizedTransmissivityFilm(double omega, std::complex<double> epsilon){

	std::complex<double> complex_exp2 (cos(2.0*BetaFunction(omega, epsilon).real())* exp(-2.0*BetaFunction(omega, epsilon).imag()), sin(2.0*BetaFunction(omega, epsilon).real())* exp(-2.0*BetaFunction(omega, epsilon).imag()));

	std::complex<double> complex_exp1 (cos(BetaFunction(omega, epsilon).real())* exp(-BetaFunction(omega, epsilon).imag()), sin(BetaFunction(omega, epsilon).real())* exp(-BetaFunction(omega, epsilon).imag()));

	return std::norm(
			(t_top_film_pPolarized(epsilon)*t_film_bottom_pPolarized(epsilon)*complex_exp1)
			/
			(1.0 + r_top_film_pPolarized(epsilon)*r_film_bottom_pPolarized(epsilon)*complex_exp2)
		);
}
 
double dielectric_function::NormalIncidenceTransmissivityFilm(double omega, std::complex<double> epsilon){

	std::complex<double> complex_exp2 (cos(2.0*BetaFunction(omega, epsilon).real())* exp(-2.0*BetaFunction(omega, epsilon).imag()), sin(2.0*BetaFunction(omega, epsilon).real())* exp(-2.0*BetaFunction(omega, epsilon).imag()));

	std::complex<double> complex_exp1 (cos(BetaFunction(omega, epsilon).real())* exp(-BetaFunction(omega, epsilon).imag()), sin(BetaFunction(omega, epsilon).real())* exp(-BetaFunction(omega, epsilon).imag()));

	return std::norm(
			(t_top_film_normalIncidence(epsilon)*t_film_bottom_normalIncidence(epsilon)*complex_exp1)
			/
			(1.0 + r_top_film_normalIncidence(epsilon)*r_film_bottom_normalIncidence(epsilon)*complex_exp2)
		);
}

double dielectric_function::NormalIncidenceTransmissivityFilm(double omega, double Te, double Tph, std::complex<double> epsilon){

	return NormalIncidenceTransmissivityFilm(omega, epsilon);
}

double dielectric_function::sPolarizedTransmissivityFilm(double omega, double Te, double Tph, std::complex<double> epsilon){
	
	return sPolarizedTransmissivityFilm(omega, epsilon);
}

double dielectric_function::pPolarizedTransmissivityFilm(double omega, double Te, double Tph, std::complex<double> epsilon){

	return pPolarizedTransmissivityFilm(omega, epsilon);
}

double dielectric_function::sPolarizedAbsorptionFilm(double omega, std::complex<double> epsilon){

	return 1.0 - (sPolarizedReflectivityFilm(omega, epsilon)+ sPolarizedTransmissivityFilm(omega, epsilon));	 
} 

double dielectric_function::pPolarizedAbsorptionFilm(double omega, std::complex<double> epsilon){

	return 1.0 - (pPolarizedReflectivityFilm(omega, epsilon)+ pPolarizedTransmissivityFilm(omega, epsilon));	 
}

double dielectric_function::NormalIncidenceAbsorptionFilm(double omega, std::complex<double> epsilon){

	//return (sPolarizedAbsorptionFilm(omega, epsilon) + pPolarizedAbsorptionFilm(omega, epsilon)) * 0.5;
	return 1.0 - (NormalIncidenceReflectivityFilm(omega, epsilon) + NormalIncidenceTransmissivityFilm(omega, epsilon)); 
}

double dielectric_function::NormalIncidenceAbsorptionFilm(double omega, double Te, double Tph, std::complex<double> epsilon){

	return  NormalIncidenceAbsorptionFilm(omega, epsilon);
}

double dielectric_function::sPolarizedAbsorptionFilm(double omega, double Te, double Tph, std::complex<double> epsilon){

	return sPolarizedAbsorptionFilm(omega, epsilon);
}

double dielectric_function::pPolarizedAbsorptionFilm(double omega, double Te, double Tph, std::complex<double> epsilon){

	return pPolarizedAbsorptionFilm(omega, epsilon);
}

/**************************************************************************
Special Properties (not Reflectivity, Transmissivity or Absorption related):
**************************************************************************/

double dielectric_function::SinglePhotonAbsorption(double omega, std::complex<double> epsilon){

	double k_ = sqrt(epsilon).imag();
	return 2.0 * (omega/hbar_eV) * k_ /speed_of_light; 	
}

double dielectric_function::SkinDepth(double omega, std::complex<double> epsilon){

	return 2.0/SinglePhotonAbsorption(omega, epsilon);
}
double dielectric_function::RealConductivity(double omega, std::complex<double> epsilon){

	return epsilon.imag() * vacuum_permittivity * (omega/hbar_eV) * 9e9; //multiply by 9e9 to convert Siemens/meter to 1/second
}

double dielectric_function::ImaginaryConductivity(double omega, std::complex<double> epsilon){
	
	return (dielectric_infty -epsilon.real()) * vacuum_permittivity * (omega/hbar_eV) * 9e9; //multiply by 9e9 to convert Siemens/meter to 1/second
}

std::complex<double> dielectric_function::ComplexConductivity(double omega, std::complex<double> epsilon){

	std::complex<double> num (RealConductivity(omega, epsilon), ImaginaryConductivity(omega, epsilon));
	return num ;
} 

// double dielectric_function::bulkReflectivity(std::complex<double> epsilon){
	
// 	std::complex<double> ref_index = std::sqrt(epsilon);
// 	double n = ref_index.real();
// 	double k = ref_index.imag();
// 	double num = std::pow((n-1),2.0) + std::pow(k, 2.0);
// 	double denom = std::pow((n+1),2.0) + std::pow(k, 2.0);
// 	return num/denom;
// }


/**************************************************************
General matrix method for n-number of films not yet implemented
**************************************************************/
/*
std::complex<double> dielectric_function::MatrixElement_m11(std::complex<double> epsilon, bool polarization){
	
	if(polarization==true){
		return std::complex<double>(0, 0);
	}
	return std::complex<double>(0, 1);
}

std::complex<double> dielectric_function::MatrixElement_m12(std::complex<double> epsilon, bool polarization){

	if(polarization==true){
		return std::complex<double>(0,0);
	}
		return std::complex<double>(0, 1);
}

std::complex<double> dielectric_function::MatrixElement_m21(std::complex<double> epsilon, bool polarization){

	if(polarization==true){
		return std::complex<double>(0, 0);
	}
		return std::complex<double>(0, 1);
}

std::complex<double> dielectric_function::MatrixElement_m22(std::complex<double> epsilon, bool polarization){

	if(polarization==true){
		return std::complex<double>(0,0);
	}
		return std::complex<double>(1,1);
}


*/



































