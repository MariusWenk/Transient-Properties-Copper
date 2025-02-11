


#ifndef FITSILAEVA_HPP
#define FITSILAEVA_HPP


#include <complex>

class fitSilaeva{


	private:
	
	double gamma0, gamma1, gamma2, gamma3, gamma4;
	double omega1, omega2, omega3, omega4;
	double f0, f1, f2, f3, f4;
	double omega_p;
	public:
		fitSilaeva();
		~fitSilaeva();

 		std::complex<double> drudeSilaeva(double);
 		std::complex<double> lorentzSilaeva1(double);
 		std::complex<double> lorentzSilaeva2(double);
 		std::complex<double> lorentzSilaeva3(double);
 		std::complex<double> lorentzSilaeva4(double);

		std::complex<double> totalSilaeva(double);
		double fitSilaevaReal(double);
		double fitSilaevaImag(double);		
		
		double _gamma0();
		double _gamma1();
		double _gamma2();
		double _gamma3();
		double _gamma4();
		
		double _omega1();
		double _omega2();
		double _omega3();
		double _omega4();
		
		double _f0();
		double _f1();
		double _f2();
		double _f3();
		double _f4();
		
		double getOmega_p();
};



#endif 
