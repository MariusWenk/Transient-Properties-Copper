#include "Laser.h"


/****************************************************
Use this constructor if you know the absorbed energy
\Delta\epsilon in J/kg converted to a.u.
*****************************************************/
Claser::Claser(){

	this ->absorbed_energy = laser_energy*(Joule_to_Hartree/kg_in_au); 

	this->pulse_duration = laser_pulse_duration * s_to_hbar_E_H;

	this->wavelength = laser_wavelength * m_to_a0;  
    	
	this->photon_energy = (planck_const_au * speed_of_light_au) / wavelength;  

    	this ->omega = photon_energy / hbar_au; 

	this ->time_peak = laser_timepeak* s_to_hbar_E_H;  // Gaussian centered at 0

	this ->material_density = material_density_const * (kg_in_au /(m_to_a0*m_to_a0*m_to_a0)); // for gold

	diel_pointer = new dielectric_function();

	boost_pointer = new BoostRootFinding1D(DOS_input);

	std::cout<<"\nThe equilibrium chemical potential mu at T = 10.000K is "<<boost_pointer->GetMuRootFinding(10000)<<" eV\n";

	drude_pointer = new dielectric_function_drude(boost_pointer);

	lorentz_pointer = new dielectric_function_lorentz(boost_pointer);
}
/****************************************************
Use this constructor if you know the incident fluence
\Phi in J/m² converted to a.u.
*****************************************************/
Claser::Claser(bool boolean){

	diel_pointer = new dielectric_function();

	boost_pointer = new BoostRootFinding1D(DOS_input);

	drude_pointer = new dielectric_function_drude(boost_pointer);

	lorentz_pointer = new dielectric_function_lorentz(boost_pointer);

	this ->boolean = boolean;  // set it to true if you know incident fluence

	this ->fluence = laser_fluence*Joule_to_Hartree/(m_to_a0*m_to_a0); 

	this->pulse_duration = laser_pulse_duration * s_to_hbar_E_H;

    	this->wavelength = laser_wavelength * m_to_a0;

    	this->photon_energy = (planck_const_au * speed_of_light_au) / wavelength;  

    	this ->omega = photon_energy / hbar_au; 

	this ->time_peak = laser_timepeak* s_to_hbar_E_H;  // Gaussian centered at 0	 	
}

 
Claser::~Claser() {

	delete diel_pointer;
	delete boost_pointer;
	delete drude_pointer;
	delete lorentz_pointer;
}

 
double Claser::GetPulseDuration() {

    	return pulse_duration;
}

double Claser::GetFluence() {

    	return fluence;
}

double Claser::GetWavelength() {

    	return wavelength;
}

double Claser::GetPhotonEnergy() {

    	return photon_energy;
}

double Claser::GetOmega() {

    	return omega;
}

double Claser::GetTimePeak(){
	
	return time_peak;
}

double Claser::GetAbsorbedEnergy(){

	return absorbed_energy;
}

double Claser::GetMaterialDensity(){

	return material_density;
}

std::complex<double> Claser::CalcDielectricFunction(double Te, double Tph){

	double energy = GetPhotonEnergy() * Hartree_to_eV; 

	double Te_conv = Te / K_in_au;

	double Tph_conv = Tph / K_in_au;

	return drude_pointer->ComplexDielectricDrude(energy, Te_conv, Tph_conv)
		+ lorentz_pointer->ComplexDielectricLorentz(energy, Te_conv);
}

std::complex<double> Claser::CalcDielectricFunction_neq(double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d){

	double energy = GetPhotonEnergy() * Hartree_to_eV; 

	double sp_density_conv = sp_density * n_au_to_SI;

	double d_density_conv = d_density * n_au_to_SI;
	
	double Tph_conv = Tph / K_in_au;
	double Te_conv = Te / K_in_au;
	double mu_sp_conv = mu_sp*Hartree_to_eV;
	double mu_d_conv = mu_d*Hartree_to_eV;

	return drude_pointer->ComplexDielectricDrude_neq(energy, Tph_conv, sp_density_conv, d_density_conv)
		+lorentz_pointer->ComplexDielectricLorentz_neq(energy, Te_conv, mu_sp_conv, mu_d_conv, d_density_conv);
}

double Claser::CalcReflectivity(double Te, double Tph){

	double energy = GetPhotonEnergy() * Hartree_to_eV; 

	return  diel_pointer->NormalIncidenceReflectivityFilm(energy, CalcDielectricFunction(Te, Tph)); 
	//return diel_pointer->pPolarizedReflectivityFilm(energy,CalcDielectricFunction(Te, Tph));
}

double Claser::CalcReflectivity_neq(double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d){

	double energy = GetPhotonEnergy() * Hartree_to_eV; 

	return diel_pointer ->NormalIncidenceReflectivityFilm(energy, CalcDielectricFunction_neq(sp_density, d_density, Te, Tph, mu_sp, mu_d));
}

double Claser::CalcTransmissivity(double Te, double Tph){

	double energy = GetPhotonEnergy() * Hartree_to_eV; 

	return  diel_pointer->NormalIncidenceTransmissivityFilm(energy, CalcDielectricFunction(Te, Tph)); 
	//return  diel_pointer->pPolarizedTransmissivityFilm(energy,CalcDielectricFunction(Te, Tph)); 
}

double Claser::CalcTransmissivity_neq(double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d){

	double energy = GetPhotonEnergy() * Hartree_to_eV; 

	return diel_pointer ->NormalIncidenceTransmissivityFilm(energy, CalcDielectricFunction_neq(sp_density, d_density, Te, Tph, mu_sp, mu_d));
}

double Claser::CalcAbsorption(double Te, double Tph){

	return 1.0 - (CalcReflectivity(Te, Tph) + CalcTransmissivity(Te, Tph));
}

double Claser::CalcAbsorption_neq(double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d){

	return 1.0 - (CalcReflectivity_neq(sp_density, d_density, Te, Tph, mu_sp, mu_d) + CalcTransmissivity_neq(sp_density, d_density, Te, Tph, mu_sp, mu_d));
}

double Claser::CalcAbsorptionCoeff(double Te, double Tph){

	double energy = GetPhotonEnergy() * Hartree_to_eV; 

	return diel_pointer->SinglePhotonAbsorption(energy,CalcDielectricFunction(Te, Tph)) / m_to_a0; // in a.u.
}

double Claser::CalcAbsorptionCoeff_neq(double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d){

	double energy = GetPhotonEnergy() * Hartree_to_eV; 

	return diel_pointer->SinglePhotonAbsorption(energy, CalcDielectricFunction_neq(sp_density, d_density, Te, Tph, mu_sp, mu_d)) / m_to_a0;
}

double Claser::LaserSourceTerm(double time){

	double S_t =  sqrt(4.0 * log(2.0) / PI)* GetAbsorbedEnergy() * GetMaterialDensity() 
		* exp(-4.0 * log(2.0) * pow(((time - GetTimePeak()) / GetPulseDuration()), 2.0)); 

	return S_t / GetPulseDuration();
}

double Claser::LaserIntensity(double time, double Te, double Tph) {

	if (boolean == true){
		double I_t =  sqrt(4.0 * log(2.0) / PI) * GetFluence() * CalcAbsorption(Te, Tph)
           	* exp(-4.0 * log(2.0) * pow(((time - GetTimePeak()) / GetPulseDuration()), 2.0)) / GetPulseDuration();

		return I_t;
	}
	std::cerr<<"\nWrong constructor or boolean\n";
	std::exit(1);
}

double Claser::LaserIntensity_neq(double time, double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d){

	if (boolean == true){

		double I_t =  sqrt(4.0 * log(2.0) / PI) * GetFluence() * CalcAbsorption_neq(sp_density, d_density, Te, Tph, mu_sp, mu_d)
           	* exp(-4.0 * log(2.0) * pow(((time - GetTimePeak()) / GetPulseDuration()), 2.0)) / GetPulseDuration();

		return I_t;
	} 
	std::cerr<<"\nWrong constructor or boolean\n";
	std::exit(1);
}


double Claser::LaserElectricField_eq(double time, double Te, double Tph) {

	double n_R = std::sqrt(CalcDielectricFunction(Te, Tph)).real();
 
	return std::sqrt(LaserIntensity(time, Te, Tph)/(speed_of_light_au*vacuum_permittivity_au*n_R));
}

double Claser::LaserElectricField_neq(double time, double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d) {

	double n_R = std::sqrt(CalcDielectricFunction_neq(sp_density, d_density, Te, Tph, mu_sp, mu_d)).real();

    	return std::sqrt( LaserIntensity_neq(time, sp_density, d_density, Te, Tph, mu_sp, mu_d)/(speed_of_light_au*vacuum_permittivity_au*n_R)); 
}

double Claser::LaserMagneticField_eq(double time, double Te, double Tph) {

    	return LaserElectricField_eq(time, Te, Tph) / speed_of_light_au;
}

double Claser::LaserMagneticField_neq(double time, double sp_density, double d_density, double Te, double Tph, double mu_sp, double mu_d) {

    	return LaserElectricField_neq(time, sp_density, d_density, Te, Tph, mu_sp, mu_d) / speed_of_light_au;
}

double Claser::KeldyshParameter(double time){

	std::cout<<"\nDo it later\n";
	return 0;
}

double Claser::KeldyshParameter(){

	std::cout<<"\nDo it later\n";
	return 0;
}


