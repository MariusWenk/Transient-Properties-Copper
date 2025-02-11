 
  
/*************************************** 
 * File:   ThreeBandWithPhonons.cpp
 * Author: Ndione, AG Rethfeld
 * 
 * Created on November 13, 2018, 10:24 AM
 ****************************************/

#include "ThreeBandWithPhonons.h"

ThreeBandWithPhonons::ThreeBandWithPhonons() : ODE(10){

	std::cerr<<"Default constructor of class \"ThreeBandWithPhonons\". Please use other constructors\n";
	std::exit(1);
}

ThreeBandWithPhonons::ThreeBandWithPhonons(std::string e_ph_filename, Cfermi_distribution *fermi_ptr) : ODE(10) {

	if (e_ph_filename != ""){

		std::ifstream data(e_ph_filename);
	
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
		std::cerr<<"\nThe file "<<e_ph_filename<< " was not correctly opened or not found! Please check the file's name or it's path\n";
		std::exit(1);
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

    	this ->laser_ptr = new Claser();
	//this ->laser_ptr = new Claser(true);

    	this ->fermi_ptr = fermi_ptr;

    	integrator_ptr = new Cintegration;

    	this ->relax_time = relaxation_time_const * s_to_hbar_E_H;

    	this ->relax_rate = 1.0 / relax_time;

	this ->inherent_linewidth = inherent_linewidth_const * eV_to_Hartree; //average

   	this ->auger_lifetime = hbar_au / inherent_linewidth;	
 
	this->auger_rate = 1.0 / auger_lifetime;	

    	this ->unit_cell_volume = unit_volume_au;
    	
	this ->fermi_energy = fermi_ptr->GetFermiEnergy();
    
	this ->start_temperature = temperature_start_const * K_in_au;
    	
	this ->n_sp_initial =  fermi_ptr->CalcEquilibriumDensity_sp(start_temperature, fermi_energy);
   
	this ->n_d_initial = fermi_ptr->CalcEquilibriumDensity_d(start_temperature, fermi_energy);

	std::cout << "\nWith the current calibration of a_sp = "<<sp_DOS_stretching<<" and a_d = "<<d_DOS_stretching
			<<" the number of sp-electrons per atom with a unit cell volume V = "<<unit_volume_SI<<" m^3 is "<<n_sp_initial*unit_volume_au
			<<" and of d-electrons is "<<n_d_initial*unit_volume_au<<" at a initial temperature of T = "<<temperature_start_const<<" K. \n";

    	this ->n_f_initial = number_of_f_electrons / unit_cell_volume;  
    	
    	this ->alpha_tot =  2.2960e7 / m_to_a0; // absorption from https://refractiveindex.info/?shelf=main&book=Au&page=Windt for \hbar\omega = 245=5.1nm eV 
    
    	this ->ph_start_temperature = temperature_start_const * K_in_au;
     
    	this ->E_sp_min = fermi_energy;
	std::cout<<"Esp = "<<E_sp_min*Hartree_to_eV<<"\n";
	
	this->d_edge = d_edge_const*eV_to_Hartree;
    
	this ->E_d_min = fermi_energy-d_edge;
	//this ->E_d_min = fermi_energy;
	std::cout<<"Ed = "<<E_d_min*Hartree_to_eV<<"\n";

   	this ->dist_sp_f = dist_sp_f_const * eV_to_Hartree;// average

	this ->E_f_min = (fermi_energy - dist_sp_f);
	std::cout<<"E_f = "<<E_f_min*Hartree_to_eV<<"\n";	
 	
	//ODEs initial conditions
    	this ->n_sp = fermi_ptr->CalcEquilibriumDensity_sp(start_temperature, fermi_energy);

    	this ->n_d = fermi_ptr->CalcEquilibriumDensity_d(start_temperature, fermi_energy);

		//std::cout<<"\nn_sp initial = "<<n_sp<<" , n_d initial = "<<n_d<<"\n";

    	this ->n_f = number_of_f_electrons/unit_volume_au;  

    	this ->n_val = fermi_ptr->CalcDensity(start_temperature, fermi_energy);

    	this ->u_val = fermi_ptr->CalcInternalEnergy(start_temperature, fermi_energy);

    	this ->T_e = start_temperature;

    	this ->T_ph = ph_start_temperature; 

    	this ->mu_sp = fermi_energy;

    	this ->mu_d = fermi_energy;

    	this ->mu_eq = fermi_energy;
    	
    	this ->t_start = simulation_t_start * s_to_hbar_E_H;
    
	this ->t_end = simulation_t_end * s_to_hbar_E_H;
    
	this ->t_step = simulation_t_step * s_to_hbar_E_H;
	
	 
	std::cout<<"\nRunning class \"ThreeBandWithPhonons\"\n";

	std::cout<<"\nSolving ODEs from t_inital = "

		<<t_start * 1e15/ s_to_hbar_E_H 

		<<" fs to t_final = "<<t_end * 1e15 / s_to_hbar_E_H 

		<<" fs with a time step of t_step = "<<t_step * 1e15 / s_to_hbar_E_H<<" fs\n";

	std::cout<<"\nVector dxdt[j] correspond to time derivatives and vector x[j] are the ODEs' solution\n";

	std::cout<<"Order of callling is important. The program is currently evaluating the following transient quantities:\n"

		<<" n_sp\n n_d\n n_f\n n_sp+d\n u_sp+d\n Te\n Tph\n mu_sp\n mu_d\n mu_eq\n";
}


ThreeBandWithPhonons::~ThreeBandWithPhonons() {

	delete [] Te_array;
	delete [] e_ph_array;
	
	delete laser_ptr;
	delete integrator_ptr;

	gsl_spline_free(e_ph_spline);
	gsl_interp_accel_free(e_ph_accel);		 
}

double ThreeBandWithPhonons::get_t_start() {

    	return t_start;
}

double ThreeBandWithPhonons::get_t_end() {
   
    	return t_end;
}

double ThreeBandWithPhonons::get_t_step() {

    	return t_step;
}

double ThreeBandWithPhonons::GetFermiEnergy() {

    	return fermi_energy;
}

double ThreeBandWithPhonons::GetStartTemperaure(){

    	return start_temperature;
}

double ThreeBandWithPhonons::get_d_edge(){
	return d_edge;
}

double ThreeBandWithPhonons::Get_n_sp_initial() {

    	return n_sp_initial;
}

double ThreeBandWithPhonons::Get_n_d_initial() {

    	return n_d_initial;
}

double ThreeBandWithPhonons::Get_n_f_initial() {

    	return n_f_initial;

}

double ThreeBandWithPhonons::Get_E_sp_min() {

    	return E_sp_min;
}

double ThreeBandWithPhonons::Get_E_d_min() {

    	return E_d_min;
}

double ThreeBandWithPhonons::Get_E_f_min() {

    	return E_f_min;
}

double ThreeBandWithPhonons::GetDistance() {

    	return dist_sp_f;
}

double ThreeBandWithPhonons::GetRelaxRate() {

    	return relax_rate;
}

double ThreeBandWithPhonons::InherentLinewidth() {

    	return inherent_linewidth;
}

double ThreeBandWithPhonons::AugerLifetime() {

    	return auger_lifetime;
}

double ThreeBandWithPhonons::GetAugerRate() {

    	return auger_rate;
}

double ThreeBandWithPhonons::GetAlphaTot() {

    	return alpha_tot;
}

double ThreeBandWithPhonons::GetUnitCellVolume() {

    	return unit_cell_volume;
}

double ThreeBandWithPhonons::GetStartingPhononTemperature(){

	return ph_start_temperature;
}

/*double ThreeBandWithPhonons::CalcAlphaTot(double Te, double Tph){

	return laser_ptr->CalcAbsorptionCoeff(Te, Tph);
}
*/
double ThreeBandWithPhonons::PhononHeatCapacity(){
  
	//constant parameter from D.R Lide, crc handbook of Chemistry and Physics (1992) 2.327MJ/Km³
	return phonon_heat_capacity_const*(Joule_to_Hartree/(K_in_au*std::pow(m_to_a0, 3.0)));
}

double ThreeBandWithPhonons::TeDependentElectronPhononCoupling(double electron_temperature){

 	if (LoadData.size() > 0){				 
		if((electron_temperature < 0.0) || (electron_temperature > Te_array[LoadData.size()-2])){
			std::cerr<<"\nThe electron temperature at which you want to interpolate is out of range. The program exits! Electron temperature: "
					<<electron_temperature<<" , T_e: "<<x[5]<<"\n";
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
	//return 4.9*((1e16*Joule_to_Hartree)/(s_to_hbar_E_H*K_in_au*std::pow(m_to_a0, 3.0)));
	//interp1D::spline interpolate;
   	//interpolate.set_points(Te_array, e_ph_array);  
	//return  interpolate(electron_temperature);
}

double ThreeBandWithPhonons::GetLossTerm(double time){
 
	return TeDependentElectronPhononCoupling(T_e) * (T_e - T_ph);		
}

double ThreeBandWithPhonons::AugerTerm(double time) {

    	return GetAugerRate() * n_d * (1.0 - (n_f / Get_n_f_initial()));
}

double ThreeBandWithPhonons::electronHoleRecombinationRate(){

	double gamma = 1e-45;
	double conversion = std::pow(m_to_a0, 6.0) / s_to_hbar_E_H;
	return gamma * conversion;
}


double ThreeBandWithPhonons::impactIonizationRate(){

	double impact_lifetime = 287e-15 * s_to_hbar_E_H;
	return 1.0 / impact_lifetime;
}

double ThreeBandWithPhonons::RelaxationTerm(double time) {

	double n_sp_eq = fermi_ptr->CalcEquilibriumDensity_sp(T_e, mu_eq);
	
    	return GetRelaxRate() * (fermi_ptr->CalcEquilibriumDensity_sp(T_e, mu_sp) - n_sp_eq);
	//return GetRelaxRate() * (n_sp - n_sp_eq);

	//return electronHoleRecombinationRate() * n_sp * (n_sp*n_sp - n_sp_eq*n_sp_eq);
	
	//double n_d_eq = fermi_ptr->CalcEquilibriumDensity_d(T_e, mu_eq);
	//return  electronHoleRecombinationRate() * n_sp * ((n_sp*n_sp)  - (n_sp_eq*n_sp_eq*n_d/n_d_eq));


	//double tem = (impactIonizationRate() * n_sp) /  (n_sp_eq * n_sp_eq);
	//return tem * (n_sp*n_sp - n_sp_eq*n_sp_eq);

}

double ThreeBandWithPhonons::AbsorbedEnergy(double time) {

	return laser_ptr->LaserSourceTerm(time);
	//return GetAlphaTot()*laser_ptr->LaserIntensity(time, T_e , T_ph);
}

double ThreeBandWithPhonons::OpticalExcitation(double time) {

	return AbsorbedEnergy(time) / laser_ptr->GetPhotonEnergy();
}

double ThreeBandWithPhonons::PotentialEnergy(double time) {

    	return (Get_E_sp_min() - Get_E_f_min()) * OpticalExcitation(time);
}

double ThreeBandWithPhonons::AugerEnergy(double time) {

    	return (Get_E_d_min() - Get_E_f_min()) * AugerTerm(time);
}

/**********************************************************
 Hands over the dynamics of conduction band electrons.
 * Optical excitation corresponds to the total absorbed 
 * energy which is divided by one photon energy. 
 * n_sp_eq correspond to the equilibrium density.
 * The relaxation time \tau_{relax} = 1.0 / relax_rate
 * is the time required to reach equilibrium density
 * 
 ***********************************************************/
void ThreeBandWithPhonons::dn_sp_dt(double time) {

    	dxdt[0] = OpticalExcitation(time) + AugerTerm(time) - RelaxationTerm(time);
	//dxdt[0] = OpticalExcitation(time) + AugerTerm(time) + impactIonization() - electronHoleRecombination();
	//int k=50;
	//dxdt[0] = OpticalExcitation(time) + (k*AugerTerm(time)) - RelaxationTerm(time);
}

/**********************************************************
 Hands over the dynamics of the valence band electrons.,
 * One electrons goes to the $4f$ shell to fill a hole, then 
 * the excess of energy may be released as a photon and extract 
 * another electron from the $d$ band. Hence we have finally lost
 * two electrons by auger process and gained one electron 
 * with the relaxation process from $sp$ to $d$ band
 **********************************************************/

void ThreeBandWithPhonons::dn_d_dt(double time) {

    	dxdt[1] = RelaxationTerm(time) - (2.0 * AugerTerm(time));
	//dxdt[1] = electronHoleRecombination()- impactIonization() - (2.0 * AugerTerm(time));
	//int k = 50;
	//dxdt[1] = RelaxationTerm(time) - ((k+1) * AugerTerm(time));
}

/**********************************************************
 Hands over the dynamics of $f$ electrons.
 * We lose electrons due to excitation by the laser
 * and gain electron from the $d$ band with three-body 
 * recombination.
 **********************************************************/
void ThreeBandWithPhonons::dn_f_dt(double time) {

    	dxdt[2] = AugerTerm(time) - OpticalExcitation(time);     
}

/***********************************************
 * Hands over the dynamics of valence electron
 * (sp+d electrons). This is needed to find 
 * the equilibrium chemical potential mu_eq
 ************************************************/

void ThreeBandWithPhonons::dn_val_dt(double time) {

    	dxdt[3] = dxdt[0] + dxdt[1];
}
/****************************************
 Hands over the total energy of valence electrons.
 We trace the kinetic energy and we add the energy due the 
 auger recombination
 ****************************************/

void ThreeBandWithPhonons::du_val_dt(double time) {

	double kinetic_energy = AbsorbedEnergy(time) - PotentialEnergy(time);

     	dxdt[4] = kinetic_energy + AugerEnergy(time)  - GetLossTerm(time);
}

/**********************************************************
 Hands over the variation of electrons temperature.
 * In this model we assume a join electrons' temperature.
 * (Interband thermalization). The temperature of the system
 * determines the equilibrum density.
 **********************************************************/

void ThreeBandWithPhonons::dTe_dt(double time) {

    	double up_faktor = (fermi_ptr->Get_CT_mu(T_e, mu_eq) * dxdt[3]) / fermi_ptr->Get_P_mu(T_e, mu_eq);
    	
	double dividend = dxdt[4] - up_faktor;

    	double down_faktor = (fermi_ptr->Get_CT_mu(T_e, mu_eq) * fermi_ptr->Get_P(T_e, mu_eq))
            	/ (fermi_ptr->Get_CT(T_e, mu_eq) * fermi_ptr->Get_P_mu(T_e, mu_eq));
    	
	double divisor = fermi_ptr->Get_CT(T_e, mu_eq) * (1.0 - down_faktor);

    	dxdt[5] = dividend / divisor;

		//std::cout<<"\n"<<dxdt[5];

}

/********************************************
*Hands over variations of phonons temperature
**********************************************/
void ThreeBandWithPhonons::dTph_dt(double time){

	dxdt[6] = (TeDependentElectronPhononCoupling(T_e) / PhononHeatCapacity()) * (T_e - T_ph);
}

/**********************************************
 * Hands over the dynamics of the chemical
 * potential of sp electrons. 
 **********************************************/
void ThreeBandWithPhonons::dmu_sp_dt(double time) {

    	double Inv_P_musp = 1.0 / fermi_ptr->Get_P_mu_sp(T_e, mu_sp);
   
	double factor_sp = fermi_ptr->Get_P_sp(T_e, mu_sp) * dxdt[5];

    	dxdt[7] = Inv_P_musp * (dxdt[0] - factor_sp);

}

/**********************************************
 * Hands over the dynamics of the chemical
 * potential of d electrons.  
 **********************************************/

void ThreeBandWithPhonons::dmu_d_dt(double time) {

    	double Inv_P_mud = 1.0 / fermi_ptr->Get_P_mu_d(T_e, mu_d);
    
	double factor_d = fermi_ptr->Get_P_d(T_e, mu_d) * dxdt[5];

    	dxdt[8] = Inv_P_mud * (dxdt[1] - factor_d);
}

/**************************************************************
 * Hands over the dynamics of the equilibrium chemical potential
 * The latter should be the same as for implicit method
 **************************************************************/

void ThreeBandWithPhonons::dmu_eq_dt(double time) {

    	double Inv_P_mu = 1.0 / fermi_ptr->Get_P_mu(T_e, mu_eq);
    
	double factor = fermi_ptr->Get_P(T_e, mu_eq) * dxdt[5];

    	dxdt[9] = Inv_P_mu * (dxdt[3] - factor);
}

/**********************************************************
*The order of calling is important in this function. 
*Please do not change unless you change your equations.
*One needs to first call the densities, 
*then the temperature and last the chemical potentials 
* **********************************************************/

void ThreeBandWithPhonons::calculate_dxdt(double time) {

	std::fill(dxdt.begin(), dxdt.end(), 0);
	
	    	dn_sp_dt(time);

		dn_d_dt(time);

    		dn_f_dt(time);

  	        dn_val_dt(time);

		du_val_dt(time);

		dTe_dt(time);

    	        dTph_dt(time);

		dmu_sp_dt(time);

	        dmu_d_dt(time);

   	        dmu_eq_dt(time);
}


