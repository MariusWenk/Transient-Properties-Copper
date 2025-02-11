#include "fitSilaeva.hpp"



fitSilaeva::fitSilaeva(){

	//fit params300 K
	/*
	this->omega_p = 13.98;
	 
	this->gamma0 = 0.06;
	this->gamma1 = 0.461;
	this->gamma2 = 1.263;
	this->gamma3 = 2.107;
	this->gamma4 = 1.493;
	
	this->omega1 = 2.53;
	this->omega2 = 3.545;
	this->omega3 = 5.261;
	this->omega4 = 7.431;
	
	this->f0 = 0.099;
	this->f1 = 0.045;
	this->f2 = 0.152;
	this->f3 = 0.146;
	this->f4 = 0.512;
	*/
	
	//fit params 10000 K
	/*
	this->omega_p = 14.26;
	this->gamma0 = 0.146;
	this->gamma1 = 0.618;
	this->gamma2 = 6.075;
	this->gamma3 = 0.072;
	this->gamma4 = 1.194;
	
	this->omega1 = 2.457;
	this->omega2 = 3.738;
	this->omega3 = 4.106;
	this->omega4 = 7.470;
	
	this->f0 = 0.184;
	this->f1 = 0.006;
	this->f2 = 0.651;
	this->f3 = 0.001;
	this->f4 = 0.368;
	*/
	
	//fit params 25000 K
	
	 this->omega_p = 15.10;
	this->gamma0 = 0.150;
	this->gamma1 = 7.133;
	this->gamma2 = 0.491;
	this->gamma3 = 0.3;
	this->gamma4 = 2.158;
	
	this->omega1 = 2.069;
	this->omega2 = 3.207;
	this->omega3 = 6.251;
	this->omega4 = 10.654;
	
	this->f0 = 0.235;
	this->f1 = 0.753;
	this->f2 = 0.013;
	this->f3 = 0.006;
	this->f4 = 0.657;
	
}

fitSilaeva::~fitSilaeva(){
}

double fitSilaeva::getOmega_p(){
	return omega_p;
}

double fitSilaeva::_f0(){
	return f0;
}

double fitSilaeva::_f1(){
	return f1;
}

double fitSilaeva::_f2(){
	return f2;
}

double fitSilaeva::_f3(){
	return f3;
}

double fitSilaeva::_f4(){
	return f4;
}

double fitSilaeva::_gamma0(){
	return gamma0;
}

double fitSilaeva::_gamma1(){
	return gamma1;
}

double fitSilaeva::_gamma2(){
	return gamma2;
}

double fitSilaeva::_gamma3(){
	return gamma3;
}

double fitSilaeva::_gamma4(){
	return gamma4;
}

double fitSilaeva::_omega1(){
	return omega1;
}

double fitSilaeva::_omega2(){
	return omega2;
}

double fitSilaeva::_omega3(){
	return omega3;
}

double fitSilaeva::_omega4(){
	return omega4;
}

std::complex<double> fitSilaeva::drudeSilaeva(double omega){

	std::complex<double> num_complex(omega*omega, omega * _gamma0());

    	return 1.0 - ((_f0()*getOmega_p()*getOmega_p()) / num_complex);
}

std::complex<double> fitSilaeva:: lorentzSilaeva1(double omega){
	double real = std::pow(_omega1(), 2.0) - std::pow(omega, 2.0);	
	
	double imag = -omega * _gamma1();
	
    	std::complex<double> num_complex(real, imag);

    	return (_f1()*getOmega_p()*getOmega_p()) / num_complex;
}

std::complex<double> fitSilaeva:: lorentzSilaeva2(double omega){
	double real = std::pow(_omega2(), 2.0) - std::pow(omega, 2.0);	
	
	double imag = -omega * _gamma2();
	
    	std::complex<double> num_complex(real, imag);

    	return (_f2()*getOmega_p()*getOmega_p()) / num_complex;
}

std::complex<double> fitSilaeva:: lorentzSilaeva3(double omega){
	double real = std::pow(_omega3(), 2.0) - std::pow(omega, 2.0);	
	
	double imag = -omega * _gamma3();
	
    	std::complex<double> num_complex(real, imag);

    	return (_f3()*getOmega_p()*getOmega_p()) / num_complex;
}

std::complex<double> fitSilaeva:: lorentzSilaeva4(double omega){
	double real = std::pow(_omega4(), 2.0) - std::pow(omega, 2.0);	
	
	double imag = -omega * _gamma4();
	
    	std::complex<double> num_complex(real, imag);

    	return (_f4()*getOmega_p()*getOmega_p()) / num_complex;
}

std::complex<double> fitSilaeva::totalSilaeva(double omega){
	return drudeSilaeva(omega) + lorentzSilaeva1(omega) + lorentzSilaeva2(omega) + lorentzSilaeva3(omega) + lorentzSilaeva4(omega);
}

double fitSilaeva::fitSilaevaReal(double omega){
	return totalSilaeva(omega).real();
}

double fitSilaeva::fitSilaevaImag(double omega){
	return totalSilaeva(omega).imag();
}









