

/*****************************************
 * File:   dielectric_function_lorentz.cpp
 * Author: Ndione, Ag Rethfeld
 *
 * Created on 27. März 2018, 15:13
 *****************************************/

#include "dielectric_function_lorentz.h"

dielectric_function_lorentz::dielectric_function_lorentz() : dielectric_function() {

	std::cerr<<"\nDefault constructor called!\nUse other constructors of class \"dielectric_function_lorentz\"\n";
	std::exit(1);
}

dielectric_function_lorentz::dielectric_function_lorentz(BoostRootFinding1D  *boost_point) : dielectric_function() {

	this ->int_pointer = new Cintegration;

    	this ->boost_point = boost_point;

	this ->unit_volume = unit_volume_SI;

	this ->fermi_energy = boost_point->GetFermiEnergy();

	this ->Te_start_temperature = temperature_start_const;

	this ->N_sp_eq = boost_point->sp_NumberofParticles(Te_start_temperature, fermi_energy);
	this ->N_d_eq = boost_point->d_NumberofParticles(Te_start_temperature, fermi_energy);
	this ->n_sp_eq = N_sp_eq / unit_volume_SI;
	this ->n_d_eq = N_d_eq / unit_volume_SI;

	this-> const_el_el = el_el_collision_const;


this -> oscillator_number = osci_num_const;

////////////////////////////////////////////////////////////////////////
//Below params obtained from a fit of 300 K experimental data of Johnson and Christy (Gold)
/*	this -> eq_strength = {4.57792525651572,4.99994587759795,12.2812487165237,24.9998286687714,34.9997412886179};
	this -> omega = {2.784,3.183,3.799,4.695,6.349};
	this -> gamma = {0.499996748382287,0.69985721035829,0.999993694756502,1.65125052470268,2.19062153152179};
*/
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//Below params obtained from a fit of 300 K experimental data of Johnson and Christy (Copper)
/*	this -> eq_strength = {4.99999895589374,4.99999723029357,14.5720346842602,24.9999972910675};
	this -> omega = {2.42094883654182,2.88348648076722,3.68041081635807,4.94292455683054};
	this -> gamma = {0.499999720648813,0.622837495786636,0.999999894705649,1.25964651529194};
*/
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//Below params obtained from a fit of 300 K experimental data of Johnson and Christy (Copper)
	this -> eq_strength = {0.923904010403028,4.99999784938873,14.9999939302167,24.9999977661911};
	this -> omega = {2.1,2.51489171228927,3.35945521942238,4.79458213409765};
	this -> gamma = {0.499999698087158,0.498769495142098,0.99999983837415,1.26596269286051};

////////////////////////////////////////////////////////////////////////

}

dielectric_function_lorentz::~dielectric_function_lorentz() {

	delete int_pointer;
}

double dielectric_function_lorentz::GetFermiEnergy(){

	return fermi_energy;
}

double dielectric_function_lorentz::Get_Te_startTemperature(){

	return Te_start_temperature;
}

double dielectric_function_lorentz::GetUnitVolume() {

    	return unit_volume;
}

double dielectric_function_lorentz::GetEqStrength1() {

    	return eq_strength[1];
}

double dielectric_function_lorentz::GetEqStrength2() {

    	return eq_strength[2];
}

double dielectric_function_lorentz::GetEqStrength3() {

    	return eq_strength[3];
}

double dielectric_function_lorentz::GetEqStrength4() {

    	return eq_strength[4];
}

double dielectric_function_lorentz::GetEqStrength5() {

    	return eq_strength[5];
}

double dielectric_function_lorentz::GetOmega1() {

    	return omega[1];
}

double dielectric_function_lorentz::GetOmega2() {

    	return omega[2];
}

double dielectric_function_lorentz::GetOmega3() {

    	return omega[3];
}

double dielectric_function_lorentz::GetOmega4() {

    	return omega[4];
}

double dielectric_function_lorentz::GetOmega5() {

    	return omega[5];
}

double dielectric_function_lorentz::GetGamma1() {

    	return gamma[1];
}

double dielectric_function_lorentz::GetGamma2() {

    	return gamma[2];
}

double dielectric_function_lorentz::GetGamma3() {

    	return gamma[3];
}

double dielectric_function_lorentz::GetGamma4() {

    	return gamma[4];
}

double dielectric_function_lorentz::GetGamma5() {

    	return gamma[5];
}

double dielectric_function_lorentz::Get_n_sp_eq() {

    	return n_sp_eq;
}

double dielectric_function_lorentz::Get_n_d_eq() {

    	return n_d_eq;
}

double dielectric_function_lorentz::Get_N_sp_eq(){

	return N_sp_eq;
}

double dielectric_function_lorentz::Get_N_d_eq(){

	return N_d_eq;
}

double dielectric_function_lorentz::get_const_el_el(){
	return const_el_el;
}

double dielectric_function_lorentz::relativeDistribution(double energy, double Te){

	double mu_relative = boost_point->GetMuRootFinding(Te)-GetFermiEnergy();

	return boost_point->GetFermiDistribution(energy, Te, mu_relative);
}

/***************************************************
We are not yet sure whether we should shift the
interband energies at \mu i.e. where the
distribution is 1/2 or whether the broadening
of the distribution should also be included
***************************************************/
double dielectric_function_lorentz::OmegaApproximation(double Te){

	return boost_point->GetMuRootFinding(Te) - (0.5*k_B_eV*Te) - (GetFermiEnergy() - 0.5*k_B_eV*Get_Te_startTemperature());
	//return (boost_point->GetMuRootFinding(Te) - GetFermiEnergy());
}

double dielectric_function_lorentz::OmegaApproximation_neq(double Te, double mu_sp){

	//return mu_sp - (0.0*k_B_eV*Te) - (boost_point->GetFermiEnergy()-(0.0*k_B_eV*Get_Te_startTemperature()));
	return mu_sp - (0.5*k_B_eV*Te) - (boost_point->GetFermiEnergy()-(0.5*k_B_eV*Get_Te_startTemperature()));
}

double dielectric_function_lorentz::dampingApproximation(double Te){

	double N_d = boost_point->d_NumberofParticles(Te, boost_point->GetMuRootFinding(Te));
	double nu_el_el = get_const_el_el() * N_d * (Get_N_d_eq() - N_d);
	return nu_el_el * hbar_eV;
}

double dielectric_function_lorentz::dampingApproximation_neq(double d_density){

	double N_d = d_density * GetUnitVolume();
	double nu_el_el = get_const_el_el() * N_d * (Get_N_d_eq() - N_d);
	return nu_el_el * hbar_eV;
}

double dielectric_function_lorentz::OmegaNonEq1(double Te){

	return GetOmega1() + OmegaApproximation(Te);
}

double dielectric_function_lorentz::OmegaNonEq1_neq(double Te, double mu_sp){

	return GetOmega1() + OmegaApproximation_neq(Te, mu_sp);
}

double dielectric_function_lorentz::OmegaNonEq2(double Te){

	return GetOmega2() + OmegaApproximation(Te);
}

double dielectric_function_lorentz::OmegaNonEq2_neq(double Te, double mu_sp){

	return GetOmega2() + OmegaApproximation_neq(Te, mu_sp);
}

double dielectric_function_lorentz::OmegaNonEq3(double Te){

	return GetOmega3() + OmegaApproximation(Te);
}

double dielectric_function_lorentz::OmegaNonEq3_neq(double Te, double mu_sp){

	return GetOmega3() + OmegaApproximation_neq(Te, mu_sp);
}

double dielectric_function_lorentz::OmegaNonEq4(double Te){

	return GetOmega4() + OmegaApproximation(Te);
}

double dielectric_function_lorentz::OmegaNonEq4_neq(double Te, double mu_sp){

	return GetOmega4() + OmegaApproximation_neq(Te, mu_sp);
}

double dielectric_function_lorentz::OmegaNonEq5(double Te){

	return GetOmega5() + OmegaApproximation(Te);
}

double dielectric_function_lorentz::OmegaNonEq5_neq(double Te, double mu_sp){

	return GetOmega5() + OmegaApproximation_neq(Te, mu_sp);
}

double dielectric_function_lorentz::GetNonEqStrength1(double Te){

	double f_Te = relativeDistribution(-GetOmega1(), Te);
	double f_Te_room = relativeDistribution(-GetOmega1(), Get_Te_startTemperature());
	return GetEqStrength1() * (f_Te / f_Te_room) ;
}

double dielectric_function_lorentz::GetNonEqStrength1_neq(double Te, double mu_d){

	double mu_relative  = mu_d-boost_point->GetFermiEnergy();
	double f_d = boost_point->GetFermiDistribution(-GetOmega1(), Te, mu_relative);
	return GetEqStrength1() * (f_d / relativeDistribution(-GetOmega1(), Get_Te_startTemperature()));
}

double dielectric_function_lorentz::GetNonEqStrength2(double Te){

	double f_Te = relativeDistribution(-GetOmega2(), Te);
	double f_Te_room = relativeDistribution(-GetOmega2(), Get_Te_startTemperature());
	return GetEqStrength2() * (f_Te / f_Te_room);
}

double dielectric_function_lorentz::GetNonEqStrength2_neq(double Te, double mu_d){

	double mu_relative = mu_d-boost_point->GetFermiEnergy();
	double f_d = boost_point->GetFermiDistribution(-GetOmega2(), Te, mu_relative);
	return GetEqStrength2() * (f_d / relativeDistribution(-GetOmega2(), Get_Te_startTemperature()));
}

double dielectric_function_lorentz::GetNonEqStrength3(double Te){

	double f_Te = relativeDistribution(-GetOmega3(), Te);
	double f_Te_room = relativeDistribution(-GetOmega3(), Get_Te_startTemperature());
	return GetEqStrength3() * (f_Te / f_Te_room);
}

double dielectric_function_lorentz::GetNonEqStrength3_neq(double Te, double mu_d){

	double mu_relative = mu_d-boost_point->GetFermiEnergy();
	double f_d = boost_point->GetFermiDistribution(-GetOmega3(), Te, mu_relative);
	return GetEqStrength3() * (f_d / relativeDistribution(-GetOmega3(), Get_Te_startTemperature()));
}

double dielectric_function_lorentz::GetNonEqStrength4(double Te){

	double f_Te = relativeDistribution(-GetOmega4(), Te);
	double f_Te_room = relativeDistribution(-GetOmega4(), Get_Te_startTemperature());
	return GetEqStrength4() * (f_Te / f_Te_room);
}

double dielectric_function_lorentz::GetNonEqStrength4_neq(double Te, double mu_d){

	double mu_relative = mu_d-boost_point->GetFermiEnergy();
	double f_d = boost_point->GetFermiDistribution(-GetOmega4(), Te, mu_relative);
	return GetEqStrength4() * (f_d / relativeDistribution(-GetOmega4(), Get_Te_startTemperature()));
}

double dielectric_function_lorentz::GetNonEqStrength5(double Te){

	double f_Te = relativeDistribution(-GetOmega5(), Te);
	double f_Te_room = relativeDistribution(-GetOmega5(), Get_Te_startTemperature());
	return GetEqStrength5() * (f_Te / f_Te_room);
}

double dielectric_function_lorentz::GetNonEqStrength5_neq(double Te, double mu_d){

	double mu_relative = mu_d-boost_point->GetFermiEnergy();
	double f_d = boost_point->GetFermiDistribution(-GetOmega5(), Te, mu_relative);
	return GetEqStrength5() * (f_d / relativeDistribution(-GetOmega5(), Get_Te_startTemperature()));
}

double dielectric_function_lorentz::GammaNonEq1(double Te){

	return GetGamma1() + dampingApproximation(Te);
}

double dielectric_function_lorentz::GammaNonEq2(double Te){

	return GetGamma2() + dampingApproximation(Te);
}

double dielectric_function_lorentz::GammaNonEq3(double Te){

	return GetGamma3() + dampingApproximation(Te);
}

double dielectric_function_lorentz::GammaNonEq4(double Te){

	return GetGamma4() + dampingApproximation(Te);
}

double dielectric_function_lorentz::GammaNonEq5(double Te){

	return GetGamma5() + dampingApproximation(Te);
}

double dielectric_function_lorentz::GammaNonEq1_neq(double d_density){

	return GetGamma1() + dampingApproximation_neq(d_density);
}

double dielectric_function_lorentz::GammaNonEq2_neq(double d_density){

	return GetGamma2() + dampingApproximation_neq(d_density);
}

double dielectric_function_lorentz::GammaNonEq3_neq(double d_density){

	return GetGamma3() + dampingApproximation_neq(d_density);
}

double dielectric_function_lorentz::GammaNonEq4_neq(double d_density){

	return GetGamma4() + dampingApproximation_neq(d_density);
}

double dielectric_function_lorentz::GammaNonEq5_neq(double d_density){

	return GetGamma5() + dampingApproximation_neq(d_density);
}

std::complex<double> dielectric_function_lorentz::ComplexLorentzOscillator_1(double omega, double Te) {

	double real = std::pow(OmegaNonEq1(Te), 2.0) - std::pow(omega, 2.0);

	double imag = -omega * GammaNonEq1(Te);

    	std::complex<double> B_complex(real, imag);

    	return GetNonEqStrength1(Te) / B_complex;
}

std::complex<double> dielectric_function_lorentz::ComplexLorentzOscillator_2(double omega, double Te) {

	double real = std::pow(OmegaNonEq2(Te), 2.0) - std::pow(omega, 2.0);

	double imag = -omega *  GammaNonEq2(Te);

    	std::complex<double> B_complex(real, imag);

    	return GetNonEqStrength2(Te) / B_complex;
}

std::complex<double> dielectric_function_lorentz::ComplexLorentzOscillator_3(double omega, double Te) {

	double real = std::pow(OmegaNonEq3(Te), 2.0) - std::pow(omega, 2.0);

	double imag = -omega * GammaNonEq3(Te);

    	std::complex<double> B_complex(real, imag);

    	return GetNonEqStrength3(Te) / B_complex;
}

std::complex<double> dielectric_function_lorentz::ComplexLorentzOscillator_4(double omega, double Te) {

	double real = std::pow(OmegaNonEq4(Te), 2.0) - std::pow(omega, 2.0);

	double imag = -omega * GammaNonEq4(Te);

    	std::complex<double> B_complex(real, imag);

    	return GetNonEqStrength4(Te) / B_complex;
}

std::complex<double> dielectric_function_lorentz::ComplexLorentzOscillator_5(double omega,double Te) {

	double real = std::pow(OmegaNonEq5(Te), 2.0) - std::pow(omega, 2.0);

	double imag = -omega * GammaNonEq5(Te);

    	std::complex<double> B_complex(real, imag);

    	return GetNonEqStrength5(Te) / B_complex;
}

std::complex<double> dielectric_function_lorentz::ComplexDielectricLorentz(double omega, double Te) {

	if(oscillator_number == 4){
		return ComplexLorentzOscillator_1(omega, Te) + ComplexLorentzOscillator_2(omega, Te)
            + ComplexLorentzOscillator_3(omega, Te) + ComplexLorentzOscillator_4(omega, Te);
	}

	if(oscillator_number == 5){
		return ComplexLorentzOscillator_1(omega, Te) + ComplexLorentzOscillator_2(omega, Te)
            + ComplexLorentzOscillator_3(omega, Te) + ComplexLorentzOscillator_4(omega, Te)
            + ComplexLorentzOscillator_5(omega, Te);
	}

	std::cerr<<"\n The Lorentz oscillator number is not possible. \n";
	return 0;
}

std::complex<double>dielectric_function_lorentz::
ComplexLorentzOscillator_1_neq(double omega, double Te, double mu_sp, double mu_d, double d_density){
	double real = std::pow(OmegaNonEq1_neq(Te, mu_sp), 2.0) - std::pow(omega, 2.0);

	double imag = -omega * GammaNonEq1_neq(d_density);

    	std::complex<double> B_complex(real, imag);

    	return GetNonEqStrength1_neq(Te, mu_d) / B_complex;
}

std::complex<double>dielectric_function_lorentz::
ComplexLorentzOscillator_2_neq(double omega, double Te, double mu_sp, double mu_d, double d_density){
	double real = std::pow(OmegaNonEq2_neq(Te, mu_sp), 2.0) - std::pow(omega, 2.0);

	double imag = -omega * GammaNonEq2_neq(d_density);

    	std::complex<double> B_complex(real, imag);

    	return GetNonEqStrength2_neq(Te, mu_d) / B_complex;
}

std::complex<double>dielectric_function_lorentz::
ComplexLorentzOscillator_3_neq(double omega, double Te, double mu_sp, double mu_d, double d_density){
	double real = std::pow(OmegaNonEq3_neq(Te, mu_sp), 2.0) - std::pow(omega, 2.0);

	double imag = -omega * GammaNonEq3_neq(d_density);

    	std::complex<double> B_complex(real, imag);

    	return GetNonEqStrength3_neq(Te, mu_d) / B_complex;
}

std::complex<double>dielectric_function_lorentz::
ComplexLorentzOscillator_4_neq(double omega, double Te, double mu_sp, double mu_d, double d_density){
	double real = std::pow(OmegaNonEq4_neq(Te, mu_sp), 2.0) - std::pow(omega, 2.0);

	double imag = -omega * GammaNonEq4_neq(d_density);

    	std::complex<double> B_complex(real, imag);

    	return GetNonEqStrength4_neq(Te, mu_d) / B_complex;
}

std::complex<double>dielectric_function_lorentz::
ComplexLorentzOscillator_5_neq(double omega, double Te, double mu_sp, double mu_d, double d_density){
	double real = std::pow(OmegaNonEq5_neq(Te, mu_sp), 2.0) - std::pow(omega, 2.0);

	double imag = -omega * GammaNonEq5_neq(d_density);

    	std::complex<double> B_complex(real, imag);

    	return GetNonEqStrength5_neq(Te, mu_d) / B_complex;
}


std::complex<double> dielectric_function_lorentz::ComplexDielectricLorentz_neq(double omega, double Te, double mu_sp, double mu_d, double d_density){

	if(oscillator_number == 4){
		return  ComplexLorentzOscillator_1_neq(omega, Te, mu_sp, mu_d, d_density)
		+ComplexLorentzOscillator_2_neq(omega, Te, mu_sp, mu_d, d_density)
		+ComplexLorentzOscillator_3_neq(omega, Te, mu_sp, mu_d, d_density)
		+ComplexLorentzOscillator_4_neq(omega, Te, mu_sp, mu_d, d_density);
	}

	if(oscillator_number == 5){
		return  ComplexLorentzOscillator_1_neq(omega, Te, mu_sp, mu_d, d_density)
		+ComplexLorentzOscillator_2_neq(omega, Te, mu_sp, mu_d, d_density)
		+ComplexLorentzOscillator_3_neq(omega, Te, mu_sp, mu_d, d_density)
		+ComplexLorentzOscillator_4_neq(omega, Te, mu_sp, mu_d, d_density)
		+ComplexLorentzOscillator_5_neq(omega, Te, mu_sp, mu_d, d_density);
	}

	std::cerr<<"\n The Lorentz oscillator number is not possible. \n";
	return 0;

}

double dielectric_function_lorentz::GetRealPart(double omega, double Te) {

    	return ComplexDielectricLorentz(omega, Te).real();
}

double dielectric_function_lorentz::GetImaginaryPart(double omega, double Te) {

    	return ComplexDielectricLorentz(omega, Te).imag();
}

double dielectric_function_lorentz::GetRealFromImaginaryPart(double omega, double Te){

	auto integrand = [ = ] (double x) ->double{

		return x * GetImaginaryPart(x, Te) / (x + omega);
	};

	double integral  = int_pointer ->integrate_qawc(integrand, 0.1, 1e10, omega);

	return 1.0 + (2.0 / PI) * integral;
}

double dielectric_function_lorentz::GetRealPart_neq(double omega, double Te, double mu_sp, double mu_d, double d_density){

	return dielectric_function_lorentz::ComplexDielectricLorentz_neq(omega, Te, mu_sp, mu_d, d_density).real();
}

double dielectric_function_lorentz::GetImaginaryPart_neq(double omega, double Te, double mu_sp, double mu_d, double d_density){

	return dielectric_function_lorentz::ComplexDielectricLorentz_neq(omega, Te, mu_sp, mu_d, d_density).imag();
}

double dielectric_function_lorentz::GetRealFromImaginaryPart_neq(double omega, double Te, double mu_sp, double mu_d, double d_density){

	auto integrand = [ = ] (double x) ->double{

		return x * GetImaginaryPart_neq(x, Te, mu_sp, mu_d, d_density) / (x + omega);
	};

	double integral  = int_pointer ->integrate_qawc(integrand, 0.1, 1e10, omega);

	return 1.0 + (2.0 / PI) * integral;
}
