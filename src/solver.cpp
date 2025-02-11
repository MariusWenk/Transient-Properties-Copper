/*********************************** 
 * File:   solver.cpp
 * Author: Ndione, Ag Rethfeld
 *
 * Created on 12. Juni 2018, 17:33
 ***********************************/  
    
#include "solver.h"


namespace pl = std::placeholders;

solver::solver(){
	
	std::cerr<<"\nDefault constructor of class \"solver\".\nPlease use other constructors\n";
	std::exit(1);
}

solver::solver(int chooseModel) {
	std::string dosFileName = DOS_input;
	this ->fermi_ptr = new Cfermi_distribution(dosFileName);

	this->laser_point = new Claser();	
	//this ->laser_point = new Claser(true);

	this ->chooseModel =  chooseModel;
 
       	Cfermi_distribution fermi_obj(dosFileName);

	std::string e_ph_filename = e_ph_coupling;
	

	std::string e_ph_filename_sp = e_ph_coupling_sp; 	//->density resolved sp-coupling only needed for 2band2T model
	std::string e_ph_filename_d = e_ph_coupling_d;  	//->density resolved  d-coupling only needed for 2band2T model

	std::cout<<"\n****************************************************************************************************"
		<<"\n\33[0;32m Used laser parameters:\n"
		<<" pulse duration = "<<laser_point->GetPulseDuration()*1e15/ s_to_hbar_E_H<<" fs FWHM\n"
		<<" photon energy = "<<laser_point->GetPhotonEnergy()*Hartree_to_eV<<" eV\n"
		<<" absorbed energy  = "<<laser_point->GetAbsorbedEnergy()*1e-6*(kg_in_au/Joule_to_Hartree)<<" MJ/kg\n"
		<<"\n\33[1;31m Make sure that you are using the right laser constructor in class Claser! \33[0m\n"
		<<"****************************************************************************************************\n";

	//if((boost::filesystem::exists(e_ph_filename_sp)) && (boost::filesystem::exists(e_ph_filename_d))){
 	if(boost::filesystem::exists(e_ph_filename)){
		if (chooseModel == 0){

			ode_point = new TwoBandWithPhonons(e_ph_filename, &fermi_obj); 
			std::cout<<"\n\33[1;31m ode_pointer pointing to class \"TwoBandWithPhonons\"\n"
				<<" Using coupling from "
				<<e_ph_filename<<" \33[0m\n";			

		}else if (chooseModel == 1){

			ode_point = new TwoBandWithPhonons2T(e_ph_filename_sp, e_ph_filename_d, &fermi_obj, 0, 1, 2);
			std::cout<<"\n\33[1;31m ode_pointer pointing to class \"TwoBandWithPhonons2T\"\n"
				<<" Using coupling from "
				<<e_ph_filename_sp<<" and "
				<<e_ph_filename_d <<" \33[0m\n";			
       		}else if (chooseModel == 2){

			ode_point = new ThreeBandWithPhonons(e_ph_filename, &fermi_obj);  
			std::cout<<"\n\33[1;31m ode_pointer pointing to class \"ThreeBandWithPhonons\"\n"
				 <<" Using coupling from "
				 <<e_ph_filename<<" \33[0m\n";

       		}else if (chooseModel == 3){

 			ode_point = new ThreeBandWithPhonons2T(e_ph_filename_sp, e_ph_filename_d, &fermi_obj, 0, 1, 2);  
			std::cout<<"\n\33[1;31m ode_pointer pointing to class \"ThreeBandWithPhonons2T\"\n"
				<<" Using coupling from "
				<<e_ph_filename_sp<<" and "
				<<e_ph_filename_d <<" \33[0m\n";			

		}else if (chooseModel == 4){

			ode_point = new FourBandWithPhonons(e_ph_filename, &fermi_obj);  
			std::cout<<"\n\33[1;31m ode_pointer pointing to class \"FourBandWithPhonons\"\n"
				 <<" Using coupling from "
				 <<e_ph_filename<<" \33[0m\n";
 
		}else if (chooseModel == 5){

			ode_point = new FourBandWithPhonons2T(e_ph_filename_sp, e_ph_filename_d, &fermi_obj, 0, 1, 2);  
			std::cout<<"\n\33[1;31m ode_pointer pointing to class \"FourBandWithPhonons2T\"\n"
				<<" Using coupling from "
				<<e_ph_filename_sp<<" and "<<e_ph_filename_d <<" \33[0m\n";			
 
		}else{
			std::cerr<<"\n\33[1;31m ode_pointer not defined. Please use index between 0 and 5 \33[0m\n";
			std::exit(1);
		}	   

	}else{
		std::cerr<<"\n\33[1;31m e_ph_Filename not found. Check the path or the file's name. The program exits. \33[0m\n";
		std::exit(1);
	}
	
	std::string appendix = file_name_appendix;
	std::string store_output_all = "var/FullDynamics"+appendix+".OUT";
	std::string store_opt_param_input = "var/FullDynamics"+appendix+".in";
	std::string store_input_params = "var/input_parameters"+appendix+".OUT";
	
 	outputfile.open(store_output_all);	//stores all variables calculated by the model
 	inputfile.open(store_opt_param_input);	//stores variables needed to calculate optical properties		
 	paramFile.open(store_input_params);  // store input data used for calculations
 	
 	if(paramFile.is_open()){
		paramFile<<"DOS file = "<<dosFileName
		<<"\ne-ph coupling = "<<e_ph_filename
		<<"\npulse duration = "<<(laser_point->GetPulseDuration()/ s_to_hbar_E_H)*1e15<<" fs FWHM"
		<<"\nmax peak of laser at "<<laser_point->GetTimePeak()/s_to_hbar_E_H*1e15<<" fs"	
		<<"\nabsorbed energy  = "<<laser_point->GetAbsorbedEnergy()*1e-6*(kg_in_au/Joule_to_Hartree)<<" MJ/kg"
		<<"\npump photon energy = "<<laser_point->GetPhotonEnergy()*Hartree_to_eV<<" eV "
		<<"Or Wavelength = "<<(laser_point->GetWavelength()/m_to_a0)*1e9<<" nm "
		<<"implies d-band absorption ratio (300K) = "<<fermi_ptr->InterbandPercentageTotal(300*K_in_au, fermi_ptr->GetFermiEnergy(), fermi_ptr->GetFermiEnergy(), laser_point->GetPhotonEnergy())
		<<"\nrelaxation time = "<<relaxation_time_const*1e15<<" fs"
		<<"\n";
	}else{
		std::cerr<<"\ninput params file "
			<<store_input_params
			<<" could not be opened. Check file path!\n Exit!\n";
		std::exit(1);
	}
	
	if(inputfile.is_open()){
		std::cout<<"\ninput for optics will be saved in "<<store_opt_param_input<<" and used if necessary!\n";
	}else{
		std::cerr<<"\n File for optics "
			<<store_opt_param_input<<" could not be opened. Check file path!\n Exit!\n";
		std::exit(1);
	}
	
	if(outputfile.is_open()){
		
		outputfile<<"#pulse duration = "<<(laser_point->GetPulseDuration()/ s_to_hbar_E_H)*1e15<<" fs FWHM\n";
		outputfile<<"#absorbed energy  = "<<laser_point->GetAbsorbedEnergy()*1e-6*(kg_in_au/Joule_to_Hartree)<<" MJ/kg\n";
		outputfile<<"#pump photon energy= "<<laser_point->GetPhotonEnergy()*Hartree_to_eV<<" eV\n";
	 	
	 	if (chooseModel == 0){
	 		outputfile<<"#time[fs]\tn_sp[1/m^3]\tn_d[1/m^3]\tn_sp+d[1/m^3]\tu[J/m^3]\tTe[K]\tTph[K]\tmu_sp[eV]\tmu_d[eV]\tmu_eq[eV]\n";
	 	}else if(chooseModel == 1){
	 		outputfile<<"\n #time[fs]\tn_sp[1/m^3]\tn_d[1/m^3]\tn_sp+d[1/m^3]\tu_sp[J/m^3]\tu_d[J/m^3]\t\tu_sp+d[J/m^3]\tTsp[K]\tTd[K]\tTeq[K]\tTph[K]\tmu_sp[eV]\tmu_d[eV]\tmu_eq[eV]\n";
	 	}else if(chooseModel == 2){
	 		outputfile<<"#time[fs]\tn_sp[1/m³]\tn_d[1/m³]\tn_f[1/m³]\tn_sp+d[1/m³]\tu[J/m³]\tTe[K]\tTph[K]\tmu_sp[eV]\tmu_d[eV]\tmu_eq[eV]\n";
	 	}else if(chooseModel == 3){
	 		outputfile<<"\n #output 3band2T\n";
	 	}else if(chooseModel == 4){
	 		outputfile<<"\n #output 4band\n";
	 	}else if(chooseModel == 5){
	 		outputfile<<"\n #output 4band2T\n";
	 	}

	}else{
		std::cerr<<"\nOutput file "
			<<store_output_all<<" could not be opened. Check file path\n!"
			<<"\nModel index also not known!" 
			<<"\nExit!\n";
		std::exit(1);
	}
	
	
	run();

	outputfile.close();
	inputfile.close();
	paramFile.close();

        std::cout << "\nClass \"solver\"ran completely\n";

}

solver::~solver() {

	delete fermi_ptr;
	delete ode_point;
	delete laser_point;
    outputfile.close();
	inputfile.close();
	paramFile.close();
}

int solver::getChosenModel(){

	return chooseModel;
}

void solver::run() {
    	
	int size = ode_point->x.size();

    	dv v(size);

    	for (int i = 0; i < size; i++) {
    
	    	v[i] = ode_point->x[i];
        }   

    	double t_start = ode_point->get_t_start();

    	double t_end = ode_point->get_t_end();

    	double t_step = ode_point->get_t_step();

    	runge_kutta_dopri5<dv> stepper1; //constant timesteps

    	controlled_runge_kutta<boost::numeric::odeint::runge_kutta_cash_karp54 < dv>> stepper2; //error-correcting time-steps

    	controlled_runge_kutta<boost::numeric::odeint::runge_kutta_dopri5 < dv >> stepper3(1e-8); //error-correcting time-steps

    	bulirsch_stoer<dv> stepper4(1e-8);

    	controlled_runge_kutta<boost::numeric::odeint::runge_kutta_fehlberg78 < dv >> stepper5(1e-7); //error-correcting time-steps

    	adams_bashforth_moulton<5, dv> stepper6; //Faster than RK dopri5 and yield same results

    	integrate_adaptive(stepper1, std::bind(&solver::odeintsolver_seb, std::ref(*this)
            ,pl::_1, pl::_2, pl::_3), v, t_start, t_end, t_step, std::bind(&solver::Print, std::ref(*this), pl::_1, pl::_2));
     
}

void solver::odeintsolver_seb(const dv &q, dv &q_dot, double time) {

	ode_point->set_new_variables(q);

    	dv changes = ode_point->get_change_vector(time);
    	
	for (unsigned int i = 0; i < q.size(); i++) {
	       
		 q_dot[i] = changes[i];
    	}
}

void solver::Print(const std::vector<double> &x, const double time) {

    	WritetoFile(x, time, getChosenModel());

	//You can print stuff calculated here for example  
 		//std::cout <<time/s_to_hbar_E_H * 1e15<<"\t"
        	//<< x[0] * n_au_to_SI<<"\t"
        	//<< x[5]/K_in_au<<"\n";
}

void solver::WritetoFile(const std::vector<double> &x, const double time, int chooseModel) {

	outputfile<< std::scientific << std::setprecision(15);
	inputfile<< std::scientific << std::setprecision(15);
 	

	if (chooseModel == 0){
		///Write to file stuff for 2-band model 1T
 	 	outputfile << (time / s_to_hbar_E_H) * 1e15 << "\t"
		<< x[0] * n_au_to_SI<<"\t"
        	<< x[1] * n_au_to_SI<<"\t"  
        	<< x[2] * n_au_to_SI<<"\t"
        	<< x[3] * Hartree_to_Joule/std::pow(a0_to_m, 3)<<"\t"
		<< x[4] / K_in_au<<"\t"
    		<< x[5] / K_in_au<<"\t"
       		<< x[6] * Hartree_to_eV << "\t"
       		<< x[7] * Hartree_to_eV << "\t"
		<< x[8] * Hartree_to_eV << "\n";

		inputfile << time / s_to_hbar_E_H * 1e15 << "\t"
		<< x[0] * n_au_to_SI<<"\t"
        	<< x[1] * n_au_to_SI<<"\t" 
        	<< x[4] /K_in_au<<"\t"       	 
    		<< x[5] / K_in_au<<"\t"
    		<< x[6] * Hartree_to_eV << "\t"
       		<< x[7] * Hartree_to_eV << "\t"
		<< x[8] * Hartree_to_eV << "\n";

	 }else if (chooseModel == 1){    		 
 		///Write to file stuff for 2-band model 2T 
 	 	outputfile << time / s_to_hbar_E_H * 1e15 << "\t"
		<< x[0] * n_au_to_SI<<"\t"
        	<< x[1] * n_au_to_SI<<"\t"  
        	<< x[2] * n_au_to_SI<<"\t"
        	<< x[3] * Hartree_to_Joule/std::pow(a0_to_m, 3)<<"\t"
		<< x[4] * Hartree_to_Joule/std::pow(a0_to_m, 3)<<"\t"
		<< x[5] * Hartree_to_Joule/std::pow(a0_to_m, 3)<<"\t"
		<< x[6] / K_in_au<<"\t"
		<< x[7] / K_in_au<<"\t"
		<< x[8] / K_in_au<<"\t"
    		<< x[9] / K_in_au<<"\t"
       		<< x[10] * Hartree_to_eV << "\t"
       		<< x[11] * Hartree_to_eV << "\t"
		<< x[12] * Hartree_to_eV << "\n";

		inputfile << time / s_to_hbar_E_H * 1e15 << "\t"
		<< x[0] * n_au_to_SI<<"\t"
        	<< x[1] * n_au_to_SI<<"\t"        	 
    		<< x[8] / K_in_au<<"\t" // check later which Te to use for optics (Tsp. Td or Teq)
    		<< x[9] / K_in_au<<"\t"	
    		<< x[10] * Hartree_to_eV << "\t"
       		<< x[11] * Hartree_to_eV << "\t"
		<< x[12] * Hartree_to_eV << "\n";	
		 
	
 	 }else if (chooseModel == 2){    		 
		///Write to file stuff for 3-band model 1T			
		outputfile <<time/s_to_hbar_E_H * 1e15<<"\t"
        	<<x[0] * n_au_to_SI<<"\t"
        	<<x[1] * n_au_to_SI<<"\t"
        	<<x[2] * n_au_to_SI<<"\t"
        	<<x[3] * n_au_to_SI<<"\t"
        	<<x[4] * Hartree_to_Joule/std::pow(a0_to_m, 3)<<"\t"
        	<<x[5] / K_in_au<<"\t"
		<<x[6] / K_in_au<<"\t" 
		<<x[7] * Hartree_to_eV<<"\t"
        	<<x[8] * Hartree_to_eV<<"\t"
        	<<x[9] * Hartree_to_eV<<"\n";
			
		inputfile <<time/s_to_hbar_E_H * 1e15<<"\t"
        	<<x[0] * n_au_to_SI<<"\t"
        	<<x[1] * n_au_to_SI<<"\t" 
        	<<x[5] / K_in_au<<"\t"    	 
		<<x[6] / K_in_au<<"\t" 
		<<x[7] * Hartree_to_eV<<"\t"
        	<<x[8] * Hartree_to_eV<<"\t"
        	<<x[9] * Hartree_to_eV<<"\n";
	
 	 }else if (chooseModel ==  3){ 
		///Write to file stuff for 3-band model 2T 
		std::cout<<"\nNot yet implemented\n";
		std::exit(1);

	 }else if(chooseModel == 4){
		///Write to file stuff for 4-band model 1T 
		outputfile <<time/s_to_hbar_E_H * 1e15<<"\t"
        	<<x[0] * n_au_to_SI<<"\t"
        	<<x[1] * n_au_to_SI<<"\t"
        	<<x[2] * n_au_to_SI<<"\t"
        	<<x[3] * n_au_to_SI<<"\t"
		<<x[4] * n_au_to_SI<<"\t"
        	<<x[5] * Hartree_to_Joule/std::pow(a0_to_m, 3)<<"\t"
        	<<x[6] / K_in_au<<"\t"
		<<x[7] / K_in_au<<"\t" 
		<<x[8] * Hartree_to_eV<<"\t"
        	<<x[9] * Hartree_to_eV<<"\t"
        	<<x[10] * Hartree_to_eV<<"\n"; 

		inputfile <<time/s_to_hbar_E_H * 1e15<<"\t"
        	<<x[0] * n_au_to_SI<<"\t"
        	<<x[1] * n_au_to_SI<<"\t"    
        	<<x[6] / K_in_au<<"\t" 	 
		<<x[7] / K_in_au<<"\t"
		<<x[8] * Hartree_to_eV<<"\t"
        	<<x[9] * Hartree_to_eV<<"\t"
        	<<x[10] * Hartree_to_eV<<"\n";
		 	 
	}else if(chooseModel == 5){
		///Write to file stuff for 4-band model 2T 
		std::cout<<"\nNot yet implemented\n";
		std::exit(1);
	 
	}else{

		std::cerr<<"\nError message. Please use an index between 0 and 5\n";
		std::exit(1);
	}
}



