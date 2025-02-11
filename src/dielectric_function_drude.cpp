
/*************************************** 
 * File:   dielectric_function_drude.h
 * Author: Ndione, Ag Rethfeld
 ***************************************
 * Created on February 9, 2017, 11:18 AM
 ***************************************/



    

#include "dielectric_function_drude.h"

dielectric_function_drude::dielectric_function_drude() : dielectric_function(){

	std::cerr<<"\nDefault constructor called!\nUse other constructors of class \"dielectric_function_drude\"\n";
	std::exit(1);
} 

dielectric_function_drude::dielectric_function_drude(BoostRootFinding1D  *boost_point) : dielectric_function(){

	this ->int_pointer = new Cintegration;   
  
    	this ->boost_point = boost_point;
 
    	this ->unit_volume = unit_volume_SI;
    	 
    	this ->fermi_energy = boost_point->GetFermiEnergy();

		this ->temperature_start = temperature_start_const;
  
    	this ->N_sp_eq = boost_point-> sp_NumberofParticles(temperature_start, fermi_energy);

    	this ->N_d_eq = boost_point-> d_NumberofParticles(temperature_start, fermi_energy);
	 
    	this ->n_sp_eq = N_sp_eq / unit_volume;

    	this ->n_d_eq = N_d_eq / unit_volume;

    	this ->sp_band_mass = sp_mass_const * m_e;

    	this ->dielectric_const_value = dielectric_infty;
    	
    	this->el_el_const = el_el_collision_const;
    	
    	this->el_ph_coll_freq_cold = el_ph_collision_const_cold;
}

dielectric_function_drude::~dielectric_function_drude() {

	delete int_pointer;
}

double dielectric_function_drude::GetBandMass_sp() {

	return sp_band_mass;
}

double dielectric_function_drude::GetUnitVolume() {

    	return unit_volume;
}

double dielectric_function_drude::GetDielectricConst() {

    	return dielectric_const_value;
}

double dielectric_function_drude::Get_n_sp_eq() {

    	return n_sp_eq;
}

double dielectric_function_drude::Get_n_d_eq() {

    	return n_d_eq;
}

double dielectric_function_drude::Get_N_d_Equilibrium() {

    	return N_d_eq;
}

double dielectric_function_drude::Get_N_sp_Equilibrium() {

    	return N_sp_eq;
}

double dielectric_function_drude::GetTemperatureStart() {

    	return temperature_start;
}

double dielectric_function_drude::GetFermiEnergy() {

    	return fermi_energy;
}

double dielectric_function_drude::get_el_el_const(){
	return el_el_const;
}

double dielectric_function_drude::get_el_ph_coll_freq_cold(){
	return el_ph_coll_freq_cold;
}

double dielectric_function_drude::GetPlasmaFreqHighTe(double Te) {

	double n_sp_Te = boost_point->sp_NumberofParticles(Te, boost_point->GetMuRootFinding(Te)) / GetUnitVolume(); 
    
	return std::sqrt((n_sp_Te * electron_charge * electron_charge) / (GetBandMass_sp() * vacuum_permittivity)) * hbar_eV;
}

double dielectric_function_drude::CalcPlasmaFreqHighTe(double sp_density) {
    
	return std::sqrt((sp_density * electron_charge * electron_charge) / (GetBandMass_sp() * vacuum_permittivity)) * hbar_eV;
}

double dielectric_function_drude::GetPlasmaFreqHighTeSquare(double Te){

	return std::pow(GetPlasmaFreqHighTe(Te), 2.0);
}

double dielectric_function_drude::CalcPlasmaFreqHighTeSquare(double sp_density){

	return std::pow(CalcPlasmaFreqHighTe(sp_density), 2.0);
}

double dielectric_function_drude::FermiVelocity(){

	return (hbar * (std::pow(3.0*M_PI*M_PI*Get_n_sp_eq(),(1.0/3.0)))) / GetBandMass_sp();
}

double dielectric_function_drude::ElectronPhononCollisionFreq(double ph_temp){

	return get_el_ph_coll_freq_cold() * (ph_temp/GetTemperatureStart()) * hbar_eV;
}

double dielectric_function_drude::ElectronElectronCollisionFreq(double Te) {
	
	double Nd0 = boost_point-> d_NumberofParticles(GetTemperatureStart(), GetFermiEnergy());

	double Nd = boost_point-> d_NumberofParticles(Te, boost_point->GetMuRootFinding(Te));

	return  get_el_el_const() * Nd * (Nd0 - Nd) * hbar_eV; // converted
}

double dielectric_function_drude::ElectronElectronCollisionFreq_neq(double d_density) {
	
	double Nd0 = boost_point-> d_NumberofParticles(GetTemperatureStart(), GetFermiEnergy());
	
	double Nd = d_density * GetUnitVolume();

	return get_el_el_const() * Nd * (Nd0 - Nd) * hbar_eV; // converted    	  
}

double dielectric_function_drude::TotalCollisionFreq(double Te, double Tph) {

	return ElectronPhononCollisionFreq(Tph) + ElectronElectronCollisionFreq(Te);
}

double dielectric_function_drude::TotalCollisionFreq_neq(double Tph, double d_density) {

	return ElectronPhononCollisionFreq(Tph) + ElectronElectronCollisionFreq_neq(d_density);
}

double dielectric_function_drude::GetTotalCollisionTime(double Te, double Tph){

	return 1.0 / TotalCollisionFreq(Te, Tph); //Be carefull with the conversion of \nu_tot 
}

std::complex<double> dielectric_function_drude::ComplexDielectricDrude(double omega, double Te, double Tph) {

    	std::complex<double> B_complex(omega*omega, omega * TotalCollisionFreq(Te, Tph));

    	return GetDielectricConst() - ((GetPlasmaFreqHighTe(Te) * GetPlasmaFreqHighTe(Te)) / B_complex);
}

std::complex<double> dielectric_function_drude::ComplexDielectricDrude_neq(double omega, double Tph, double sp_density, double d_density) {

    	std::complex<double> B_complex(omega*omega, omega * TotalCollisionFreq_neq(Tph, d_density));

    	return GetDielectricConst() - ((CalcPlasmaFreqHighTe(sp_density) * CalcPlasmaFreqHighTe(sp_density)) / B_complex);
}	

std::complex<double> dielectric_function_drude::ComplexDielectricDrude(double omega, double Te, double Tph, double wavevector) {

    	return ComplexDielectricDrude(omega, Te, Tph);
}

double dielectric_function_drude::GetRealPart(double omega, double Te, double Tph) {

    	return ComplexDielectricDrude(omega, Te, Tph).real();
}

double dielectric_function_drude::GetRealPart_neq(double omega, double Tph, double sp_density, double d_density) {

    	return ComplexDielectricDrude_neq(omega, Tph, sp_density, d_density).real();
}

double dielectric_function_drude::GetRealPart(double omega, double Te, double Tph, double wavevector) {

    	return GetRealPart(omega, Te, Tph);
}

double dielectric_function_drude::GetImaginaryPart(double omega, double Te, double Tph) {

    	return ComplexDielectricDrude(omega, Te, Tph).imag();
}

double dielectric_function_drude::GetImaginaryPart_neq(double omega, double Tph, double sp_density, double d_density) {

    	return ComplexDielectricDrude_neq(omega, Tph, sp_density, d_density).imag();
}

double dielectric_function_drude::GetImaginaryPart(double omega, double Te, double Tph, double wavevector) {

    	return GetImaginaryPart(omega, Te, Tph);
}

double dielectric_function_drude::GetRealFromImaginaryPart(double omega, double Te, double Tph){

	auto integrand = [ = ] (double x) ->double {

		return x * GetImaginaryPart(x, Te, Tph) / (x + omega);
	};

	double integral = int_pointer -> integrate_qawc(integrand, 0.1, 1e10, omega);

	return 1.0 + (2.0 / PI) * integral;
}

double dielectric_function_drude::GetRealFromImaginaryPart_neq(double omega, double Tph, double sp_density, double d_density){

	auto integrand = [ = ] (double x) ->double {

		return x * GetImaginaryPart_neq(x, Tph, sp_density, d_density) / (x + omega);
	};

	double integral = int_pointer -> integrate_qawc(integrand, 0.1, 1e10, omega);

	return 1.0 + (2.0 / PI) * integral;
}

double dielectric_function_drude::GetRealFromImaginaryPart(double omega, double Te, double Tph, double wavevector){

	 return GetRealFromImaginaryPart(omega, Te, Tph);
}

double dielectric_function_drude::dcConductivity(double Te, double Tph){

	double sigma_0 = vacuum_permittivity * (GetPlasmaFreqHighTeSquare(Te)/hbar_eV/hbar_eV) / (TotalCollisionFreq(Te, Tph)/hbar_eV) ;

	return sigma_0;// * 9e9; //multiply by 9e9 to convert S/m in Hz
}

/***************************
Conductivity related output:
***************************/

double dielectric_function_drude::dcConductivity_neq(double Tph, double sp_density, double d_density){

	double sigma_0 = vacuum_permittivity * (CalcPlasmaFreqHighTeSquare(sp_density)/hbar_eV/hbar_eV) / (TotalCollisionFreq_neq(Tph, d_density)/hbar_eV);
	return sigma_0;// * 9e9;  //multiply by 9e9 to convert S/m in Hz
}

double dielectric_function_drude::Resistivity(double Te, double Tph){

	return 1.0 / dcConductivity(Te, Tph);
}

double dielectric_function_drude::Resistivity_neq(double Tph, double sp_density, double d_density){
	
	return 1.0  / dcConductivity_neq(Tph, sp_density, d_density);
}












