 
/**************************************** 
 * File:   FourBandWithPhonons.cpp	*
 * Author: Ndione, AG Rethfeld		*
 *					*
 * Created on May 06, 2020, 01:08 PM
 ***************************************/

#include "FourBandWithPhonons.h"

FourBandWithPhonons::FourBandWithPhonons() : ODE(11){

	std::cerr<<"Default constructor of class \"FourBandWithPhonons\". Please use other constructors\n";
	std::exit(1);
}

FourBandWithPhonons::FourBandWithPhonons(std::string filename, Cfermi_distribution *fermi_ptr) : ODE(11) {

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
		}
	
	}else{
		std::cerr<<"\nThe file "<<filename<< " was not correctly opened or not found! Please check the file's name or it's path\n";
		std::exit(1);
	}

    	this->Te_array = new double [LoadData.size()];

	this->e_ph_array = new double [LoadData.size()];

	for (unsigned int i = 0; i != LoadData.size(); i++){

		Te_array[i] = LoadData[i][0]; //*(K_in_au);//Convert Te in a.u. ///convert only if SI coupling is used as input
		
		e_ph_array[i] = LoadData[i][1]; //*((Joule_to_Hartree)/(s_to_hbar_E_H*K_in_au*std::pow(m_to_a0, 3.0)));//Convert G in a.u. ///convert only if SI coupling is used as input

	} 
	 
	e_ph_accel = gsl_interp_accel_alloc();

	e_ph_spline = gsl_spline_alloc(gsl_interp_akima, (int) LoadData.size());

	gsl_spline_init(e_ph_spline, Te_array, e_ph_array, (int) LoadData.size());

    	this ->laser_ptr = new Claser();
	//this ->laser_ptr = new Claser(true);

    	this ->fermi_ptr = fermi_ptr;

    	integrator_ptr = new Cintegration;

    	this ->relax_time = relaxation_time_const * s_to_hbar_E_H;

    	this ->relax_rate = 1.0 / relax_time;

	this ->inherent_linewidth = inherent_linewidth_const * eV_to_Hartree; //average

   	this ->lifetime_4f_states = hbar_au / inherent_linewidth;	
 
	this->rate_4f_states = 1.0 / lifetime_4f_states;	

    	this ->unit_cell_volume = unit_volume_au;
    	
	this ->fermi_energy = fermi_ptr->GetFermiEnergy();
    
	this ->start_temperature = temperature_start_const * K_in_au;
    	
	this ->n_6sp_initial =  fermi_ptr->CalcEquilibriumDensity_sp(start_temperature, fermi_energy);
   
	this ->n_5d_initial = fermi_ptr->CalcEquilibriumDensity_d(start_temperature, fermi_energy);

	std::cout << "\nWith the current calibration of a_sp = "<<sp_DOS_stretching<<" and a_d = "<<d_DOS_stretching
			<<" the number of sp-electrons per atom with a unit cell volume V = "<<unit_volume_SI<<" m^3 is "<<n_6sp_initial*unit_volume_au
			<<" and of d-electrons is "<<n_5d_initial*unit_volume_au<<" at a initial temperature of T = "<<temperature_start_const<<" K. \n";

	this ->n_5p_initial = number_of_p_electrons / unit_cell_volume;

    	this ->n_4f_initial = number_of_f_electrons / unit_cell_volume;  
    	
    	this ->alpha_tot = 5.28985e7 / m_to_a0; //for \hbar\omega = 85 eV=14.58638 nm // not really needed here

    	this ->alpha_inter = alpha_tot;
    
    	this ->ph_start_temperature = temperature_start_const * K_in_au;
     
    	this ->E_6sp_min = fermi_energy;
    
	this ->E_5d_min = fermi_energy;

	this ->dist_6sp_5p = dist_sp_p_const * eV_to_Hartree;// average splitted bands

	this ->E_5p_min = fermi_energy - dist_6sp_5p;

   	this ->dist_6sp_4f = dist_sp_f_const * eV_to_Hartree;// average splitted bands

	this ->E_4f_min = fermi_energy - dist_6sp_4f;	
 
	this ->lifetime_5p_states = lifetime_p_states_const * s_to_hbar_E_H;//= 0.078e-15 * s_to_hbar_E_H; Guess

	this ->rate_5p_states = 1.0 / lifetime_5p_states;

	//ODEs initial conditions
	this ->n_6sp = fermi_ptr->CalcEquilibriumDensity_sp(start_temperature, fermi_energy);

    	this ->n_5d = fermi_ptr->CalcEquilibriumDensity_d(start_temperature, fermi_energy);

    	this ->n_5p = n_5p_initial;  

    	this ->n_4f = n_4f_initial;  

    	this ->n_val = fermi_ptr->CalcDensity(start_temperature, fermi_energy);

    	this ->u_val = fermi_ptr->CalcInternalEnergy(start_temperature, fermi_energy);

    	this ->T_e = start_temperature;

    	this ->T_ph = ph_start_temperature; 

    	this ->mu_6sp = fermi_energy;

    	this ->mu_5d = fermi_energy;

    	this ->mu_eq = fermi_energy;

    	this ->t_start = simulation_t_start * s_to_hbar_E_H;
    
	this ->t_end = simulation_t_end * s_to_hbar_E_H;
    
	this ->t_step = simulation_t_step * s_to_hbar_E_H;
	
	 
	std::cout<<"\nRunning class \"FourBandWithPhonons\"\n";

	std::cout<<"\nSolving ODEs from t_inital = "

		<<t_start * 1e15/ s_to_hbar_E_H 

		<<" fs to t_final = "<<t_end * 1e15 / s_to_hbar_E_H 

		<<" fs with a time step of t_step = "<<t_step * 1e15 / s_to_hbar_E_H<<"fs\n";

	std::cout<<"\nVector dxdt[j] correspond to time derivatives and vector x[j] are the ODEs' solution\n";

	std::cout<<"Order of callling is important. The program is currently evaluating the following transient quantities:\n"

		<<" n_6sp\n n_5d\n n_5p\n n_4f\n n_sp+d\n u_sp+d\n Te\n Tph\n mu_6sp\n mu_5d\n mu_eq\n";

}

FourBandWithPhonons::~FourBandWithPhonons() {

	delete [] Te_array;
	delete [] e_ph_array;
	
	delete laser_ptr;
	delete integrator_ptr;

	gsl_spline_free(e_ph_spline);
	gsl_interp_accel_free(e_ph_accel);		 
}

double FourBandWithPhonons::get_t_start() {

    	return t_start;
}

double FourBandWithPhonons::get_t_end() {
   
    	return t_end;
}

double FourBandWithPhonons::get_t_step() {

    	return t_step;
}

double FourBandWithPhonons::GetFermiEnergy() {

    	return fermi_energy;
}

double FourBandWithPhonons::GetStartTemperaure(){

    	return start_temperature;
}

double FourBandWithPhonons::Get_n_6sp_initial() {

    	return n_6sp_initial;
}

double FourBandWithPhonons::Get_n_5d_initial() {

    	return n_5d_initial;
}

double FourBandWithPhonons::Get_n_4f_initial() {

    	return n_4f_initial;
}

double FourBandWithPhonons::Get_E_6sp_min() {

    	return E_6sp_min;
}

double FourBandWithPhonons::Get_E_5d_min() {

    	return E_5d_min;
}

double FourBandWithPhonons::Get_E_5p_min() {

    	return E_5p_min;
}

double FourBandWithPhonons::Get_E_4f_min() {

    	return E_4f_min;
}

double FourBandWithPhonons::GetDistance6sp_4f() {

    	return dist_6sp_4f;
}

double FourBandWithPhonons::GetDistance6sp_5p() {

    	return dist_6sp_5p;
}

double FourBandWithPhonons::GetRelaxRate() {

    	return relax_rate;
}

double FourBandWithPhonons::InherentLinewidth() {

    	return inherent_linewidth;
}

double FourBandWithPhonons::GetLifetime4fStates() {

    	return lifetime_4f_states;
}

double FourBandWithPhonons::GetRate4fStates() {

    	return rate_4f_states;
}

double FourBandWithPhonons::GetAlphaTot() {

    	return alpha_tot;
}

double FourBandWithPhonons::GetAlphaInter() {

    	return alpha_inter;
}

double FourBandWithPhonons::GetUnitCellVolume() {

    	return unit_cell_volume;
}

double FourBandWithPhonons::GetStartingPhononTemperature(){

	return ph_start_temperature;

}

double FourBandWithPhonons::Get_n_5p_iniial(){

	return n_5p_initial;
}

double FourBandWithPhonons::GetLifetime5pStates(){

	return lifetime_5p_states;
}

double FourBandWithPhonons::GetRate5pStates(){

	return rate_5p_states;
}   

double FourBandWithPhonons::PhononHeatCapacity(){
  
	//constant parameter from D.R Lide, crc handbook of Chemistry and Physics (1992) 2.327MJ/Km³
	return phonon_heat_capacity_const*(Joule_to_Hartree/(K_in_au*std::pow(m_to_a0, 3.0)));
}

double FourBandWithPhonons::TeDependentElectronPhononCoupling(double electron_temperature){

 	if (LoadData.size() > 0){		
		 
		if((electron_temperature < 0.0) || (electron_temperature > Te_array[LoadData.size()-2])){

			std::cerr<<"\nThe electron temperature at which you want to interpolate is out of range. The program exits!\n";
			std::exit(1);

		}else{
			if((electron_temperature > 0.0) && (electron_temperature < Te_array[0])){

				return e_ph_array[0];
				 
			}else{		
				return gsl_spline_eval(e_ph_spline, electron_temperature, e_ph_accel);
			}			
		}
	
	}else {	
        	std::cerr << "\nNo points found in the data you want to interpolate. The program exits!\n";
        	std::exit(1);
	}  
}

double FourBandWithPhonons::GetLossTerm(){
 
	return TeDependentElectronPhononCoupling(T_e) * (T_e - T_ph);		
}

double FourBandWithPhonons::AugerTerm_5pBand() {

    	return GetRate5pStates() * n_5d * (1.0 - (n_5p / Get_n_5p_iniial()));
}

double FourBandWithPhonons::AugerTerm_4fBand() {

    	return GetRate4fStates() * n_5p * (1.0 - (n_4f / Get_n_4f_initial()));
}

double FourBandWithPhonons::RelaxationTerm() {
 
    	return GetRelaxRate() * (fermi_ptr->CalcEquilibriumDensity_sp(T_e, mu_6sp) - fermi_ptr->CalcEquilibriumDensity_sp(T_e, mu_eq));
}

double FourBandWithPhonons::AbsorbedEnergy(double time) {

	return laser_ptr->LaserSourceTerm(time);
}

double FourBandWithPhonons::OpticalExcitation(double time) {

	return AbsorbedEnergy(time) / laser_ptr->GetPhotonEnergy();
}

double FourBandWithPhonons::PotentialEnergy(double time) {

    	return (Get_E_6sp_min() - Get_E_4f_min()) * OpticalExcitation(time);
}

double FourBandWithPhonons::AugerEnergy5p_4f() {

    	return (Get_E_5p_min() - Get_E_4f_min()) * AugerTerm_4fBand();
}

double FourBandWithPhonons::AugerEnergy5d_5p() {

    	return (Get_E_5d_min() - Get_E_5p_min()) * AugerTerm_5pBand();
}

void FourBandWithPhonons::dn_6sp_dt(double time) {

    	dxdt[0] = OpticalExcitation(time) + AugerTerm_4fBand() + AugerTerm_5pBand() - RelaxationTerm();
}

void FourBandWithPhonons::dn_5d_dt(double time) {

    	dxdt[1] = RelaxationTerm() - (2.0 * AugerTerm_5pBand());
}

void FourBandWithPhonons::dn_5p_dt(double time) {

    	dxdt[2] = AugerTerm_5pBand() - (2.0 * AugerTerm_4fBand());
}

void FourBandWithPhonons::dn_4f_dt(double time) {

    	dxdt[3] = AugerTerm_4fBand() - OpticalExcitation(time);
}

void FourBandWithPhonons::dn_val_dt(double time) {

    	dxdt[4] = dxdt[0] + dxdt[1];
}

void FourBandWithPhonons::du_val_dt(double time) {

	double kinetic_energy = AbsorbedEnergy(time) - PotentialEnergy(time) ;
	
	double Auger_energy5p4f = (Get_E_5p_min() - Get_E_4f_min()) * AugerTerm_4fBand(); 

	double Auger_energy5d5p = (Get_E_5d_min() - Get_E_5p_min()) * AugerTerm_5pBand();

	double coupling_e_ph = GetLossTerm();

    	dxdt[5] = kinetic_energy + Auger_energy5p4f + Auger_energy5d5p - coupling_e_ph;
}

void FourBandWithPhonons::dTe_dt(double time) {

    	double up_faktor = (fermi_ptr->Get_CT_mu(T_e, mu_eq) * dxdt[4]) / fermi_ptr->Get_P_mu(T_e, mu_eq);
    	
	double dividend = dxdt[5] - up_faktor;

    	double down_faktor = (fermi_ptr->Get_CT_mu(T_e, mu_eq) * fermi_ptr->Get_P(T_e, mu_eq))
            	/ (fermi_ptr->Get_CT(T_e, mu_eq) * fermi_ptr->Get_P_mu(T_e, mu_eq));
    	
	double divisor = fermi_ptr->Get_CT(T_e, mu_eq) * (1.0 - down_faktor);

    	dxdt[6] = dividend / divisor;
}

void FourBandWithPhonons::dTph_dt(double time) {

    	dxdt[7] = (TeDependentElectronPhononCoupling(T_e) / PhononHeatCapacity()) * (T_e - T_ph);
}

void FourBandWithPhonons::dmu_6sp_dt(double time) {

    	double Inv_P_musp = 1.0 / fermi_ptr->Get_P_mu_sp(T_e, mu_6sp);
   
	double factor_sp = fermi_ptr->Get_P_sp(T_e, mu_6sp) * dxdt[6];

    	dxdt[8] = Inv_P_musp * (dxdt[0] - factor_sp);
}

void FourBandWithPhonons::dmu_5d_dt(double time) {

    	double Inv_P_mud = 1.0 / fermi_ptr->Get_P_mu_d(T_e, mu_5d);
    
	double factor_d = fermi_ptr->Get_P_d(T_e, mu_5d) * dxdt[6];

    	dxdt[9] = Inv_P_mud * (dxdt[1] - factor_d);
}

void FourBandWithPhonons::dmu_eq_dt(double time) {

    	double Inv_P_mu = 1.0 / fermi_ptr->Get_P_mu(T_e, mu_eq);
    
	double factor = fermi_ptr->Get_P(T_e, mu_eq) * dxdt[6];

    	dxdt[10] = Inv_P_mu * (dxdt[4] - factor);
}

void FourBandWithPhonons::calculate_dxdt(double time) {

	std::fill(dxdt.begin(), dxdt.end(), 0);
	
	    	dn_6sp_dt(time);

		dn_5d_dt(time);

		dn_5p_dt(time);

    		dn_4f_dt(time);

  	        dn_val_dt(time);

		du_val_dt(time);

		dTe_dt(time);

    	        dTph_dt(time);

		dmu_6sp_dt(time);

	        dmu_5d_dt(time);

   	        dmu_eq_dt(time);
	 
}








































