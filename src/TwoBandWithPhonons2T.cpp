/* ************************************
 * File:   TwoBandWithPhonons2T.hpp
 * Author: ndione
 *
 * Created on June 09, 2020, 3:21 PM
 *************************************/

#include "TwoBandWithPhonons2T.hpp"

TwoBandWithPhonons2T::TwoBandWithPhonons2T() : ODE(13) {

	std::cerr<<"Default constructor of class \"TwoBandWithPhonons2T\". Please use other constructors\n";
	std::exit(1);
}

TwoBandWithPhonons2T::TwoBandWithPhonons2T(std::string sp_filename2D, std::string d_filename2D, Cfermi_distribution *fermi_ptr,  uint col0, uint col1, uint col2) : ODE(13) {

	this ->col0 = col0;

	this ->col1 = col1;

	this ->col2 = col2;

	sp_interp_pointer2D = new ufd::Interpolation2D(sp_filename2D, col0, col1, col2);

	d_interp_pointer2D = new ufd::Interpolation2D(d_filename2D, col0, col1, col2);

	//laser_ptr = new Claser(true);
	laser_ptr = new Claser();

    	this ->fermi_ptr = fermi_ptr;

    	integrator_ptr = new Cintegration;

    	this ->elec_start_temperature = temperature_start_const*K_in_au;

    	this ->ph_start_temperature = temperature_start_const*K_in_au;

    	this ->fermi_energy = fermi_ptr->GetFermiEnergy();

    	this ->relax_time = relaxation_time_const * s_to_hbar_E_H;

    	this ->relax_rate = 1.0 / relax_time;

    	this ->unit_cell_volume = unit_volume_au;   
      
    	this ->n_sp_initial = fermi_ptr->CalcEquilibriumDensity_sp(elec_start_temperature, fermi_energy);

    	this ->n_d_initial = fermi_ptr->CalcEquilibriumDensity_d(elec_start_temperature, fermi_energy);

		std::cout << "\nWith the current calibration of a_sp = "<<sp_DOS_stretching<<" and a_d = "<<d_DOS_stretching
			<<" the number of sp-electrons per atom with a unit cell volume V = "<<unit_volume_SI<<" m^3 is "<<n_sp_initial*unit_volume_au
			<<" and of d-electrons is "<<n_d_initial*unit_volume_au<<" at a initial temperature of T = "<<temperature_start_const<<" K. \n";

    	this ->total_density_2_bands = fermi_ptr->CalcDensity(elec_start_temperature, fermi_energy);

	this ->interband_time_2T = interband_time_2T_const*s_to_hbar_E_H;
  	
	//ODE's initial conditions
    	this ->n_sp = fermi_ptr->CalcEquilibriumDensity_sp(elec_start_temperature, fermi_energy);

    	this ->n_d = fermi_ptr->CalcEquilibriumDensity_d(elec_start_temperature, fermi_energy); 

    	this ->n_val = fermi_ptr->CalcDensity(elec_start_temperature, fermi_energy);  
	
	this ->u_sp = fermi_ptr->CalcInternalEnergy_sp(elec_start_temperature, fermi_energy);

	this ->u_d = fermi_ptr->CalcInternalEnergy_d(elec_start_temperature, fermi_energy);

    	this ->u_val = fermi_ptr->CalcInternalEnergy(elec_start_temperature, fermi_energy);

    	this ->T_sp = elec_start_temperature;

	this ->T_d = elec_start_temperature;

	this ->T_eq = elec_start_temperature;

    	this ->T_ph = ph_start_temperature;

    	this ->mu_sp = fermi_energy;                                

    	this ->mu_d = fermi_energy;

    	this ->mu_eq = fermi_energy;

    	this ->t_start = simulation_t_start * s_to_hbar_E_H;  

    	this ->t_end = simulation_t_end * s_to_hbar_E_H;

    	this ->t_step = simulation_t_step * s_to_hbar_E_H;

	std::cout<<"\nRunning class \"TwoBandWithPhonons2T\"\n";

	std::cout<<"\nSolving ODEs from t_inital = "

		<<t_start * 1e15/ s_to_hbar_E_H 

		<<" fs to t_final = "<<t_end * 1e15 / s_to_hbar_E_H 

		<<" fs with a time step of t_step = "<<t_step * 1e15 / s_to_hbar_E_H<<"fs\n";

	std::cout<<"\nVector dxdt[j] correspond to time derivatives and vector x[j] are the ODEs' solution\n";

	std::cout<<"Order of callling is important. The program is currently evaluating the following transient quantities:\n"

		<<" n_sp\n n_d\n n_sp+d\n u_sp\n u_d\n u_sp+d\n T_sp\n T_d\n T_eq\n Tph\n mu_sp\n mu_d\n mu_eq\n";
	
}

TwoBandWithPhonons2T::~TwoBandWithPhonons2T() {
	
	delete laser_ptr;
	delete integrator_ptr;
	delete sp_interp_pointer2D;	
	delete d_interp_pointer2D;	
}

double TwoBandWithPhonons2T::get_t_start() {

    	return t_start;
}

double TwoBandWithPhonons2T::get_t_end() {

    	return t_end;
}

double TwoBandWithPhonons2T::get_t_step() {

    	return t_step;
}

double TwoBandWithPhonons2T::GetStartTemperature(){

    	return elec_start_temperature;
}

double TwoBandWithPhonons2T::GetStartTemperature_ph(){

    	return ph_start_temperature;
}
 
double TwoBandWithPhonons2T::GetFermiEnergy() {

    	return fermi_energy;
}

double TwoBandWithPhonons2T::Get_n_sp_initial() {

    	return n_sp_initial;
}
 
double TwoBandWithPhonons2T::Get_n_d_initial() {

    	return n_d_initial;
}

double TwoBandWithPhonons2T::GetRelaxationTime() {

    	return relax_time;
}

double TwoBandWithPhonons2T::GetRelaxRate() {

    	return relax_rate;
}

double TwoBandWithPhonons2T::GetUnitCellVolume() {

    	return unit_cell_volume;
}

double TwoBandWithPhonons2T::GetTotalDensity() {

    	return total_density_2_bands;
}

double TwoBandWithPhonons2T::CalcAlphaTot(double Te, double Tph){

	return laser_ptr->CalcAbsorptionCoeff(Te, Tph);
}

double TwoBandWithPhonons2T::GetInterbandTime2T(){

	return interband_time_2T;
}

double TwoBandWithPhonons2T::PhononHeatCapacity(){
  
	//constant parameter from D.R Lide, crc handbook of Chemistry and Physics (1992) 2.327MJ/Km³
	return phonon_heat_capacity_const*(Joule_to_Hartree/(K_in_au*std::pow(m_to_a0, 3.0)));
	//oder Cl = 2.5MJ/Km³A. M. James and M. P. Lord, Maxmillian’s Chemical and Physical Data (Maxmillian Press, London, 1972)
}

double TwoBandWithPhonons2T::TeDensityResolvedCoupling_sp(double sp_num_electron, double Te){

	return sp_interp_pointer2D->getValue(sp_num_electron, Te);
}

double TwoBandWithPhonons2T::TeDensityResolvedCoupling_d(double d_num_electron, double Te){

	return d_interp_pointer2D->getValue(d_num_electron, Te);
}

double TwoBandWithPhonons2T::spCoupling(){

	double N_sp = n_sp * GetUnitCellVolume();

	//return  TeDensityResolvedCoupling_sp(N_sp, T_sp);
	return sp_interp_pointer2D->getValue(N_sp, T_sp);
}

double TwoBandWithPhonons2T::dCoupling(){	
	
	double N_d = n_d * unit_volume_au; 

	//interpolation data starts with 10 d-electron sharp. Hence we need to round our number of electrons.
	// We cut it after 3 digits for numerical convinience
	double N_d_rounded = floorf(N_d * 1e3) / 1e3;

	//return  TeDensityResolvedCoupling_d(N_d_rounded, T_d);
	return d_interp_pointer2D->getValue(N_d_rounded, T_d);
}

double TwoBandWithPhonons2T::TotalCoupling(){

	return spCoupling() + dCoupling();
}

double TwoBandWithPhonons2T::GetLossTerm_sp(){

	return spCoupling() * (T_sp - T_ph);	 	 	
}

double TwoBandWithPhonons2T::GetLossTerm_d(){

	return dCoupling() * (T_d - T_ph);	 	 	
}

double TwoBandWithPhonons2T::GetLossTerm_tot(){

	return GetLossTerm_sp() + GetLossTerm_d();	 	 	
}

double TwoBandWithPhonons2T::ElectronElectronCoupling_sp(){

	return fermi_ptr->Get_CT_sp(T_sp, mu_sp) / GetInterbandTime2T();
}

double TwoBandWithPhonons2T::ElectronElectronCoupling_d(){

	return fermi_ptr->Get_CT_d(T_d, mu_d) / GetInterbandTime2T();
}

double TwoBandWithPhonons2T::ElectronElectronCoupling_tot(){

	return ElectronElectronCoupling_sp() +  ElectronElectronCoupling_d(); 
	/*double num = fermi_ptr->Get_Ce_sp(T_sp, mu_sp) * fermi_ptr->Get_Ce_d(T_d, mu_d);
	double denom = fermi_ptr->Get_Ce_sp(T_sp, mu_sp) + fermi_ptr->Get_Ce_d(T_d, mu_d);
	return (num / denom) * (1.0/GetInterbandTime2T());*/
}

double TwoBandWithPhonons2T::InterbandPercentageFromDynamics() {

    	return fermi_ptr->NumberStatesTotal_d(T_d, mu_sp, mu_d, laser_ptr->GetPhotonEnergy())
            / (fermi_ptr->NumberStatesTotal_sp(T_sp, mu_sp, mu_d, laser_ptr->GetPhotonEnergy()) + fermi_ptr->NumberStatesTotal_d(T_d, mu_sp, mu_d, laser_ptr->GetPhotonEnergy()));
}

double TwoBandWithPhonons2T::IntrabandPercentageFromDynamics() {

    	return fermi_ptr->NumberStatesTotal_sp(T_sp, mu_sp, mu_d, laser_ptr->GetPhotonEnergy())
            / (fermi_ptr->NumberStatesTotal_sp(T_sp, mu_sp, mu_d, laser_ptr->GetPhotonEnergy()) + fermi_ptr->NumberStatesTotal_d(T_d, mu_sp, mu_d, laser_ptr->GetPhotonEnergy()));
	//return 1.0 - InterbandPercentageFromDynamics();
}


double TwoBandWithPhonons2T::OpticalExcitation(double time) {
 	
	/*return (InterbandPercentageFromDynamics() 
		* laser_ptr->CalcAbsorptionCoeff_neq(T_ph, n_sp, n_d) 
		* laser_ptr->LaserIntensity_neq(time, T_ph, n_sp, n_d)) 
		/ laser_ptr->GetPhotonEnergy();
	*/
	return (InterbandPercentageFromDynamics() *  laser_ptr->LaserSourceTerm(time)) / laser_ptr->GetPhotonEnergy();
}

double TwoBandWithPhonons2T::RelaxationTerm() {

    	return GetRelaxRate() * (fermi_ptr->CalcEquilibriumDensity_sp(T_sp, mu_sp) - fermi_ptr->CalcEquilibriumDensity_sp(T_eq, mu_eq));
}

double TwoBandWithPhonons2T::AbsorbedEnergy(double time) {

	//return laser_ptr->CalcAbsorptionCoeff_neq(T_ph, n_sp, n_d) *  laser_ptr->LaserIntensity_neq(time, T_ph, n_sp, n_d);
	return laser_ptr->LaserSourceTerm(time);
}

void TwoBandWithPhonons2T::dn_sp_dt(double time) {

    	dxdt[0] = OpticalExcitation(time) - RelaxationTerm();
}

void TwoBandWithPhonons2T::dn_d_dt(double time) {

    	dxdt[1] = -dxdt[0];
}

void TwoBandWithPhonons2T::dn_val_dt(double time) {

    	dxdt[2] = dxdt[0] + dxdt[1];
}

void TwoBandWithPhonons2T::du_sp_dt(double time) {

	//double S_intra = IntrabandPercentageFromDynamics() * AbsorbedEnergy(time);
	double S_intra =  AbsorbedEnergy(time);

	double coupling_ee = ElectronElectronCoupling_tot() * (T_d - T_sp);

	double coupling_eph = GetLossTerm_sp(); 

	dxdt[3] =  S_intra + coupling_ee - coupling_eph;
}

void TwoBandWithPhonons2T::du_d_dt(double time) {

	//double S_inter = InterbandPercentageFromDynamics() * AbsorbedEnergy(time);
	double S_inter = 0.0;

	double coupling_ee = ElectronElectronCoupling_tot() * (T_d - T_sp);

	double coupling_eph = GetLossTerm_d(); 

	dxdt[4] =  S_inter - coupling_ee - coupling_eph;
}

void TwoBandWithPhonons2T::du_val_dt(double time) {

	dxdt[5] = dxdt[3] + dxdt[4];
}

void TwoBandWithPhonons2T::dT_sp_dt(double time) {

	double num_factor = (1.0 / fermi_ptr->Get_P_mu_sp(T_sp, mu_sp))
            * fermi_ptr->Get_CT_mu_sp(T_sp, mu_sp)
            * dxdt[0];

    	double dividend = dxdt[3] - num_factor;

    	double div_factor = (fermi_ptr->Get_CT_mu_sp(T_sp, mu_sp) * fermi_ptr->Get_P_sp(T_sp, mu_sp))
            / (fermi_ptr->Get_CT_sp(T_sp, mu_sp) * fermi_ptr->Get_P_mu_sp(T_sp, mu_sp));

    	double divisor = fermi_ptr->Get_CT_sp(T_sp, mu_sp) * (1.0 - div_factor);

    	dxdt[6] = dividend / divisor;   	
}

void TwoBandWithPhonons2T::dT_d_dt(double time) {

    	double num_factor = (1.0 / fermi_ptr->Get_P_mu_d(T_d, mu_d))
            * fermi_ptr->Get_CT_mu_d(T_d, mu_d)
            * dxdt[1];

    	double dividend = dxdt[4] - num_factor;

    	double div_factor = (fermi_ptr->Get_CT_mu_d(T_d, mu_d) * fermi_ptr->Get_P_d(T_d, mu_d))
            / (fermi_ptr->Get_CT_d(T_d, mu_d) * fermi_ptr->Get_P_mu_d(T_d, mu_d));

    	double divisor = fermi_ptr->Get_CT_d(T_d, mu_d) * (1.0 - div_factor);

    	dxdt[7] = dividend / divisor;	
}

void TwoBandWithPhonons2T::dT_eq_dt(double time) {
	
	double num_factor = (1.0 / fermi_ptr->Get_P_mu(T_eq, mu_eq))
            * fermi_ptr->Get_CT_mu(T_eq, mu_eq)
            * dxdt[2];

    	double dividend = dxdt[5] - num_factor;

    	double div_factor = (fermi_ptr->Get_CT_mu(T_eq, mu_eq) * fermi_ptr->Get_P(T_eq, mu_eq))
            / (fermi_ptr->Get_CT(T_eq, mu_eq) * fermi_ptr->Get_P_mu(T_eq, mu_eq));

    	double divisor = fermi_ptr->Get_CT(T_eq, mu_eq) * (1.0 - div_factor);

    	dxdt[8] = dividend / divisor;
}

void TwoBandWithPhonons2T::dTph_dt(double time){

	dxdt[9] = GetLossTerm_tot() / PhononHeatCapacity();
}

void TwoBandWithPhonons2T::dmu_sp_dt(double time) {

	double Inv_P_musp = 1.0 / fermi_ptr->Get_P_mu_sp(T_sp, mu_sp);

    	double factor_sp = fermi_ptr->Get_P_sp(T_sp, mu_sp) * dxdt[6];

    	dxdt[10] = Inv_P_musp * (dxdt [0] - factor_sp);
}

void TwoBandWithPhonons2T::dmu_d_dt(double time) {

	double Inv_P_mud = 1.0 / fermi_ptr->Get_P_mu_d(T_d, mu_d);

    	double factor_d = fermi_ptr->Get_P_d(T_d, mu_d) * dxdt[7];

    	dxdt[11] = Inv_P_mud * (dxdt[1] - factor_d);
  
}

void TwoBandWithPhonons2T::dmu_eq_dt(double time) {

	double Inv_P_mu = 1.0 / fermi_ptr->Get_P_mu(T_eq, mu_eq);

    	double factor = fermi_ptr->Get_P(T_eq, mu_eq) * dxdt[8];

    	dxdt[12] = Inv_P_mu * (dxdt[2] - factor);
    	
}

void TwoBandWithPhonons2T::calculate_dxdt(double time) {

    	std::fill(dxdt.begin(), dxdt.end(), 0); 

    	dn_sp_dt(time);

    	dn_d_dt(time);

	dn_val_dt(time);

	du_sp_dt(time);

	du_d_dt(time);

    	du_val_dt(time);

    	dT_sp_dt(time);

	dT_d_dt(time);

	dT_eq_dt(time);

    	dTph_dt(time);

    	dmu_sp_dt(time);

    	dmu_d_dt(time);

    	dmu_eq_dt(time);
}











































