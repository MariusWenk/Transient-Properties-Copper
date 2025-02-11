/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TwoBandWithPhonons.cpp
 * Author: ndione
 * 
 * Created on August 21, 2018, 3:17 PM
 */

#include "TwoBandWithPhonons.h"

TwoBandWithPhonons::TwoBandWithPhonons() : ODE(9) {

	std::cerr<<"Default constructor of class \"TwoBandWithPhonons\". Please use other constructors\n";
	std::exit(1);
}

TwoBandWithPhonons::TwoBandWithPhonons(std::string filename,  Cfermi_distribution *fermi_ptr) : ODE(9) {

	if (filename != ""){

		std::ifstream data(filename);
	
		std::vector <double> value(2);

		if (data.is_open()){

			while(data >>value[0] >> value[1]){
			
				if(LoadData.size()>0){
		
					value[0] *= 1;    //convert later
	
					value[1] *= 1;

					LoadData.push_back(value);

				}else{

					LoadData.push_back(value);

				}
			
			}
		}else{
			std::cerr<<"\nThe file "<<filename<< " was not correctly opened or not found! Please check the file's name or it's path\n";
			std::exit(1);
		}  
	
	}

        this->Te_array = new double [LoadData.size()];

	this->e_ph_array = new double [LoadData.size()];

	for (unsigned int i = 0; i != LoadData.size(); i++){

		Te_array[i] = LoadData[i][0]; //*(K_in_au);//Convert Te in a.u. ///convert only if SI coupling is used as input
		
		e_ph_array[i] = LoadData[i][1]; //*((Joule_to_Hartree)/(s_to_hbar_E_H*K_in_au*std::pow(m_to_a0, 3.0)));//Convert G in a.u. ///convert only if SI coupling is used as input
		
	} 
		 
	//if((sizeof(Te_array)/(sizeof(*Te_array))) != (sizeof(e_ph_array)/sizeof(*e_ph_array))){
	// 
	//		std::cerr<<"\nThe data you want to interpolate do not have the same sizes. The program exits!\n";
	//		std::exit(1);
	//} 

	e_ph_accel = gsl_interp_accel_alloc();

	e_ph_spline = gsl_spline_alloc(gsl_interp_akima, (int) LoadData.size());

	gsl_spline_init(e_ph_spline, Te_array, e_ph_array, (int) LoadData.size());

	//interp_pointer1D = new ufd::Interpolation1D(filename, 0, 1);
    
    	laser_ptr = new Claser();
	//laser_ptr = new Claser(true);

    	this ->fermi_ptr = fermi_ptr;

    	integrator_ptr = new Cintegration;

    	this ->elec_start_temperature = temperature_start_const*K_in_au;

    	this ->ph_start_temperature = temperature_start_const*K_in_au;

    	this ->fermi_energy = fermi_ptr->GetFermiEnergy();

    	this ->relax_time = relaxation_time_const*s_to_hbar_E_H;

    	this ->relax_rate = 1.0 / relax_time;

    	this ->unit_cell_volume = unit_volume_au;   
      
    	this ->n_sp_initial = fermi_ptr->CalcEquilibriumDensity_sp(elec_start_temperature, fermi_energy);

    	this ->n_d_initial = fermi_ptr->CalcEquilibriumDensity_d(elec_start_temperature, fermi_energy);

		std::cout << "\nWith the current calibration of a_sp = "<<sp_DOS_stretching<<" and a_d = "<<d_DOS_stretching
			<<" the number of sp-electrons per atom with a unit cell volume V = "<<unit_volume_SI<<" m^3 is "<<n_sp_initial*unit_volume_au
			<<" and of d-electrons is "<<n_d_initial*unit_volume_au<<" at a initial temperature of T = "<<temperature_start_const<<" K. \n";
		
    	this ->total_density_2_bands = fermi_ptr->CalcDensity(elec_start_temperature, fermi_energy);
   
	 
    	//this ->alpha_tot = 7.7089e7 / m_to_a0;        // total absorption from https://refractiveindex.info/?shelf=main&book=Au&page=Johnson \hbar \omega = 1.55 eV  
    	this ->alpha_tot = 6.1355e7 / m_to_a0; // total absorption from https://refractiveindex.info/?shelf=main&book=Au&page=Johnson \hbar \omega = 3.1 eV  	

	this -> twoPhotonAbsorptionCoeff = 0.1e-5 * ((m_to_a0 * s_to_hbar_E_H)/Joule_to_Hartree);// please convert me

	//ODEs initial conditions
    	this ->n_sp = fermi_ptr->CalcEquilibriumDensity_sp(elec_start_temperature, fermi_energy);

    	this ->n_d = fermi_ptr->CalcEquilibriumDensity_d(elec_start_temperature, fermi_energy); 

    	this ->n_val = fermi_ptr->CalcDensity(elec_start_temperature, fermi_energy);  

    	this ->u_val = fermi_ptr->CalcInternalEnergy(elec_start_temperature, fermi_energy);

    	this ->T_e = elec_start_temperature;

    	this ->T_ph = ph_start_temperature;

    	this ->mu_sp = fermi_energy;

    	this ->mu_d = fermi_energy;

    	this ->mu_eq = fermi_energy;

    	this ->t_start = simulation_t_start * s_to_hbar_E_H;  //t_sart could be -3*FWHM. Here we start at -5000fs just for the convolution

    	this ->t_end = simulation_t_end * s_to_hbar_E_H;

    	this ->t_step = simulation_t_step * s_to_hbar_E_H;

	std::cout<<"\nRunning class \"TwoBandWithPhonons\"\n";

	std::cout<<"\nSolving ODEs from t_inital = "

		<<t_start * 1e15 / s_to_hbar_E_H 

		<<" fs to t_final = "<<t_end * 1e15 / s_to_hbar_E_H 

		<<" fs with a time step of t_step = "<<t_step * 1e15 / s_to_hbar_E_H<<" fs\n";

	std::cout<<"\nVector dxdt[j] correspond to time derivatives and vector x[j] are the ODEs' solution\n";

	std::cout<<"Order of callling is important. The program is currently evaluating the following transient quantities:\n"

		<<" n_sp\n n_d\n n_sp+d\n u_sp+d\n Te\n Tph\n mu_sp\n mu_d\n mu_eq\n";

}


TwoBandWithPhonons::~TwoBandWithPhonons() {
	
	delete laser_ptr;
	delete integrator_ptr;
	
	delete Te_array;
	delete e_ph_array;

	gsl_spline_free(e_ph_spline);
	gsl_interp_accel_free(e_ph_accel);	
}

double TwoBandWithPhonons::get_t_start() {

    	return t_start;
}

double TwoBandWithPhonons::get_t_end() {

    	return t_end;
}

double TwoBandWithPhonons::get_t_step() {

    	return t_step;
}

double TwoBandWithPhonons::GetStartTemperature(){

    	return elec_start_temperature;
}

double TwoBandWithPhonons::GetStartTemperature_ph(){

    	return ph_start_temperature;
}
 
double TwoBandWithPhonons::GetFermiEnergy() {

    	return fermi_energy;
}

double TwoBandWithPhonons::Get_n_sp_initial() {

    	return n_sp_initial;
}
 
double TwoBandWithPhonons::Get_n_d_initial() {

    	return n_d_initial;
}

double TwoBandWithPhonons::GetRelaxationTime() {

    	return relax_time;
}

double TwoBandWithPhonons::GetRelaxRate() {

    	return relax_rate;
}

double TwoBandWithPhonons::GetAlphaTot() {

    	return alpha_tot;
}
 
double TwoBandWithPhonons::GetUnitCellVolume() {

    	return unit_cell_volume;
}

 
double TwoBandWithPhonons::GetTotalDensity() {

    	return total_density_2_bands;
}

double TwoBandWithPhonons::GetTwoPhotonAbsorptionCoeff(){

	return twoPhotonAbsorptionCoeff;
}

double TwoBandWithPhonons::PhononHeatCapacity(){
  
	//constant parameter from D.R Lide, crc handbook of Chemistry and Physics (1992) 2.327MJ/Km³
	return phonon_heat_capacity_const*(Joule_to_Hartree/(K_in_au*std::pow(m_to_a0, 3.0)));
}

double TwoBandWithPhonons::TeDependentElectronPhononCoupling(double electron_temperature){

 	if (LoadData.size() > 0){	
	
		if((electron_temperature < Te_array[0]) || (electron_temperature > Te_array[LoadData.size()-2])){
	
			std::cerr<<"\nThe electron temperature at which you want to interpolate is out of range. The program exits!\n";
			std::exit(1);
		}else{
			
			return gsl_spline_eval(e_ph_spline, electron_temperature, e_ph_accel);			 
			//return 0;

		}
	
	}else {
	
        	std::cerr << "\nNo points found in the data you want to interpolate. The program exits!\n";
        	std::exit(1);

	}	
}


double TwoBandWithPhonons::GetLossTerm(){
 
	return TeDependentElectronPhononCoupling(T_e) * (T_e - T_ph);
}

double TwoBandWithPhonons::InterbandPercentageFromDynamics() {

	return fermi_ptr->InterbandPercentageTotal(T_e, mu_sp, mu_d, laser_ptr->GetPhotonEnergy());
}

double TwoBandWithPhonons::IntrabandPercentageFromDynamics() {

	return fermi_ptr->IntrabandPercentageTotal(T_e, mu_sp, mu_d, laser_ptr->GetPhotonEnergy());
}

double TwoBandWithPhonons::electronHoleRecombinationRate(){

	double gamma = 1e-44;
	double conversion = std::pow(m_to_a0, 6.0) / s_to_hbar_E_H;
	return gamma * conversion;
}

double TwoBandWithPhonons::impactIonizationRate(){

	double impact_lifetime = 28.6e-15 * s_to_hbar_E_H;
	return 1.0 / impact_lifetime;
}

double TwoBandWithPhonons::RelaxationTerm() {

	double n_sp_eq = fermi_ptr->CalcEquilibriumDensity_sp(T_e, mu_eq);

    	return GetRelaxRate() * (fermi_ptr->CalcEquilibriumDensity_sp(T_e, mu_sp) - n_sp_eq);
	//return GetRelaxRate() * (n_sp - n_sp_eq);

	//return electronHoleRecombinationRate() * n_sp * (n_sp*n_sp - n_sp_eq*n_sp_eq);
	
	//double n_d_eq = fermi_ptr->CalcEquilibriumDensity_d(T_e, mu_eq);
	//double n_sp_eq = fermi_ptr->CalcEquilibriumDensity_sp(T_e, mu_eq);
	//return  electronHoleRecombinationRate() * n_sp * ((n_sp*n_sp)  - (n_sp_eq*n_sp_eq*n_d/n_d_eq));

	
	//double tem = (impactIonizationRate() * n_sp) /  (n_sp_eq * n_sp_eq);
	//return tem * (n_sp*n_sp - n_sp_eq*n_sp_eq);
}

double TwoBandWithPhonons::AbsorbedEnergy(double time) {

	return laser_ptr->LaserSourceTerm(time);

	//return laser_ptr->CalcAbsorptionCoeff_neq(n_sp, n_d, T_e, T_ph, mu_eq) *  laser_ptr->LaserIntensity_neq(time, n_sp, n_d, T_e, T_ph, mu_eq);
	//at the moment not yet clear whether to use mu_eq, mu_sp or mu_d for noneq optics
	//return laser_ptr->CalcAbsorptionCoeff(T_e, T_ph) *  laser_ptr->LaserIntensity(time, T_e, T_ph);
}

double TwoBandWithPhonons::OpticalExcitation(double time) {
 	
	return (InterbandPercentageFromDynamics() *  AbsorbedEnergy(time)) / laser_ptr->GetPhotonEnergy();
}
double TwoBandWithPhonons::TwoPhotonAbsorption(double time){

	//double I_t_square = -->implement later;

	//return  GetTwoPhotonAbsorptionCoeff() * I_t_square; 
	std::cerr<<"\nNot yet well tested.\nExits!\n";
	std::exit(1);
	return 0.0;
}

double TwoBandWithPhonons::TwoPhotonAbsorptionFunction(double time){

	return  TwoPhotonAbsorption(time) / (2.0 * laser_ptr->GetPhotonEnergy());
}

/**********************************************************
 *Hands over the dynamics of conduction band electrons.
 * Electrons are excited from d to sp-band
 ***********************************************************/

void TwoBandWithPhonons::dn_sp_dt(double time) {

    	dxdt[0] = OpticalExcitation(time) - RelaxationTerm();
	//dxdt[0] = OpticalExcitation(time) + TwoPhotonAbsorptionFunction(time)  - RelaxationTerm(time);
	//dxdt[0] = OpticalExcitation(time) + impactIonization() - electronHoleRecombination();
}

/**********************************************************
 *Hands over the dynamics of the valence band electrons.,
 * Electrons relax from the sp to the d-band
 ***********************************************************/
void TwoBandWithPhonons::dn_d_dt(double time) {

    	dxdt[1] = -dxdt[0];
}

/********************************************************
 *  Hands over the total electrons density of conduction
 *and valence band, which we denote here dn_valence
 ********************************************************/
void TwoBandWithPhonons::dn_val_dt(double time) {

    	dxdt[2] = dxdt[0] + dxdt[1];
}

/****************************************************
 *Hands over the valence electrons' internal energy. For the
 *two band, this correspond to the total energy we put into the system 
 ******************************************************/
void TwoBandWithPhonons::du_val_dt(double time) {

    	dxdt[3] = AbsorbedEnergy(time) - GetLossTerm();
	//dxdt[3] = AbsorbedEnergy(time) + TwoPhotonAbsorption(time) - GetLossTerm(time);
}

/**********************************************************
 Hands over the dynamics of electrons temperature.
 *Since we assume that the energy exange is fast, 
 *hence, sp and d-bands have the same temperature 
 **********************************************************/
void TwoBandWithPhonons::dTe_dt(double time) {

    	double num_factor = (1.0 / fermi_ptr->Get_P_mu(T_e, mu_eq))
            * fermi_ptr->Get_CT_mu(T_e, mu_eq)
            * dxdt[2];

    	double dividend = dxdt[3] - num_factor;

    	double div_factor = (fermi_ptr->Get_CT_mu(T_e, mu_eq) * fermi_ptr->Get_P(T_e, mu_eq))
            / (fermi_ptr->Get_CT(T_e, mu_eq) * fermi_ptr->Get_P_mu(T_e, mu_eq));

    	double divisor = fermi_ptr->Get_CT(T_e, mu_eq) * (1.0 - div_factor);

    	dxdt[4] = dividend / divisor;
}

/********************************************
*Hands over variations of phonons temperature
**********************************************/
void TwoBandWithPhonons::dTph_dt(double time){

   	dxdt[5] = (TeDependentElectronPhononCoupling(T_e) / PhononHeatCapacity()) * (T_e - T_ph);
}


/***********************************************************
 *Hands over the dynamics of sp electrons' chemical potential
 ***********************************************************/
void TwoBandWithPhonons::dmu_sp_dt(double time) {

    	double Inv_P_musp = 1.0 / fermi_ptr->Get_P_mu_sp(T_e, mu_sp);

    	double factor_sp = fermi_ptr->Get_P_sp(T_e, mu_sp) * dxdt[4];

    	dxdt[6] = Inv_P_musp * (dxdt [0] - factor_sp);
}

/***********************************************************
 *Hands over the dynamics of d electrons' chemical potential
 ***********************************************************/
void TwoBandWithPhonons::dmu_d_dt(double time) {

    	double Inv_P_mud = 1.0 / fermi_ptr->Get_P_mu_d(T_e, mu_d);

    	double factor_d = fermi_ptr->Get_P_d(T_e, mu_d) * dxdt[4];

    	dxdt[7] = Inv_P_mud * (dxdt[1] - factor_d);

		//if(!std::isfinite(dxdt[7]) and fermi_ptr->Get_P_mu_d(T_e, mu_d) == 0 and std::abs(dxdt[1] - factor_d) <= 1e-20){
		//	dxdt[7] = 0;
		//}
}

/***********************************************************
 *Hands over the dynamics of equilibrium chemical potential
 ***********************************************************/
void TwoBandWithPhonons::dmu_eq_dt(double time) {

    	double Inv_P_mu = 1.0 / fermi_ptr->Get_P_mu(T_e, mu_eq);

    	double factor = fermi_ptr->Get_P(T_e, mu_eq) * dxdt[4];

    	dxdt[8] = Inv_P_mu * (dxdt[2] - factor);
}

/********************************************************************
The order of calling is  important in this function. 
*One needs to start with the densities n_sp, n_d, n_val and then u_val. 
*Since the temperature is needed for the chemical potential, one should
*then call Te, mu_sp, mu_d and mu_eq
**********************************************************************/
void TwoBandWithPhonons::calculate_dxdt(double time) {

    	std::fill(dxdt.begin(), dxdt.end(), 0); 

    	dn_sp_dt(time);

    	dn_d_dt(time);

    	dn_val_dt(time);

    	du_val_dt(time);

    	dTe_dt(time);

    	dTph_dt(time);

    	dmu_sp_dt(time);

    	dmu_d_dt(time);

    	dmu_eq_dt(time);
}










