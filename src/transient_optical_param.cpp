/****************************************
 * File:   transient_optical_param.cpp
 * Author: Ndione, Ag Rethfeld
 *
 * Created on 12. Juni 2018, 17:33
*****************************************/

#include "transient_optical_param.h"

transient_optical_param::transient_optical_param() {
	std::cerr<<"\nDefault constructor called!\nUse other constructors of class transient_optical_param\n";
	std::exit(1);
}

transient_optical_param::transient_optical_param(std::string input_optics_file,
	dielectric_function_drude *drude_point, dielectric_function_lorentz *lorentz_point)
{

	this->drude_point = drude_point;

	this->lorentz_point = lorentz_point;

	this ->wavelength = probe_wavelength_const; // probe wavelength

    	this ->omega_probe = 2.0*M_PI*hbar_eV*speed_of_light/wavelength;  // probe frequency in eV

    	std::vector<double> values(8);

    	std::ifstream data(input_optics_file);

    	if (data.is_open()) {

		while (data >> values[0] >> values[1] >> values[2] >> values[3] >> values[4] >> values[5] >> values[6] >> values[7]) {

			opt_param.push_back(values);
        	}

    	}else{

        std::cerr << "\nThe file " << input_optics_file<<" that you have called was not found or not correctly opened." << "\n"
	          <<" Please check the filename or it's path :)" << "\n";
        std::exit(1);
    	}

    	this ->n_size =  opt_param.size();

    	this ->Time = new double [(int) n_size];

    	this ->n_sp = new double [(int) n_size];

    	this ->n_d = new double [(int) n_size];

	this ->T_electron = new double [(int) n_size];

	this ->T_lattice = new double [(int) n_size];

    	this->mu_sp = new double [(int)  n_size];

    	this->mu_d = new double [(int)  n_size];

	this ->mu_eq = new double [(int) n_size];

    	for (unsigned int j = 0; j <  n_size; j++) {

        	Time[j] = opt_param[j][0];
        	n_sp[j] = opt_param[j][1];
        	n_d[j] = opt_param[j][2];
		T_electron[j] = opt_param[j][3];
		T_lattice[j] = opt_param[j][4];
		mu_sp[j] = opt_param[j][5];
		mu_d[j] = opt_param[j][6];
		mu_eq[j] = opt_param[j][7];

    	}

	this ->complexDielectricFunction = new std::complex<double> [n_size];
	this ->n_complex= new std::complex<double> [n_size];
	this ->sigma_r= new double [n_size];
	this ->sigma_i= new double [n_size];
	this ->dcConductivity= new double [n_size];
	this ->e_ph_nu= new double [n_size];
    	this ->spd_nu= new double [n_size];
    	this ->tot_nu= new double [n_size];

    	this ->R_n_bulk= new double [n_size];
    	this ->R_s_polarized_bulk= new double [n_size];
    	this ->R_p_polarized_bulk= new double [n_size];

	this ->T_n_bulk= new double [n_size];
	this ->T_s_polarized_bulk= new double [n_size];
	this ->T_p_polarized_bulk= new double [n_size];

	this ->A_n_bulk= new double [n_size];
	this ->A_s_polarized_bulk= new double [n_size];
	this ->A_p_polarized_bulk= new double [n_size];

	this ->R_n_film= new double [n_size];
    	this ->R_s_polarized_film= new double [n_size];
    	this ->R_p_polarized_film= new double [n_size];

    	this ->T_n_film= new double [n_size];
	this ->T_s_polarized_film= new double [n_size];
	this ->T_p_polarized_film= new double [n_size];

	this ->A_n_film= new double [n_size];
	this ->A_p_polarized_film= new double [n_size];
	this ->A_s_polarized_film= new double [n_size];

	this -> probe_wavelength_amount = probe_wavelength_amount_const;
	this -> min_probe_wavelength = min_probe_wavelength_const;
	this -> max_probe_wavelength = max_probe_wavelength_const;

	this -> time_step_factor = probe_time_step_factor;

	this->complexDielectricFunction_eq = new std::complex<double> [probe_wavelength_amount];
	this->omega = new double [probe_wavelength_amount];
	this->T_film_eq = new double [probe_wavelength_amount];
	this->R_film_eq = new double [probe_wavelength_amount];

}

transient_optical_param::~transient_optical_param() {

    	delete [] Time;
    	delete [] n_sp;
    	delete [] n_d;
	delete [] T_lattice;
	delete [] T_electron;
	delete [] mu_sp;
	delete [] mu_d;
	delete [] mu_eq;

    	delete [] complexDielectricFunction;
    	delete [] n_complex;
	delete [] sigma_r;
	delete [] sigma_i;
	delete [] dcConductivity;
    	delete [] e_ph_nu;
    	delete [] spd_nu;
    	delete [] tot_nu;
    	delete [] R_n_bulk;
    	delete [] R_s_polarized_bulk;
    	delete [] R_p_polarized_bulk;
	delete [] T_n_bulk;
	delete [] T_s_polarized_bulk;
	delete [] T_p_polarized_bulk;
	delete [] A_n_bulk;
	delete [] A_s_polarized_bulk;
	delete [] A_p_polarized_bulk;
	delete [] R_n_film;
    	delete [] R_s_polarized_film;
    	delete [] R_p_polarized_film;
    	delete [] T_n_film;
	delete [] T_s_polarized_film;
	delete [] T_p_polarized_film;
	delete [] A_n_film;
	delete [] A_p_polarized_film;
	delete [] A_s_polarized_film;
}

double transient_optical_param::GetOmegaProbe() {

    	return omega_probe;
}

double transient_optical_param::GetWavelength() {

    	return wavelength;
}

void transient_optical_param::OpticalPropertiesDrudeLorentz() {

	std::string appendix = file_name_appendix;


	Cfermi_distribution *fermi_ptr = new Cfermi_distribution(DOS_input);


	/*
	std::string transientOptics = "var/optics/TransientOptics"+appendix+"_"+std::to_string((int)(probe_wavelength_const*1e9))+".OUT";
	std::ofstream optics(transientOptics);
	if (optics.is_open()) {
		optics<<std::scientific<<std::setprecision(15);
		optics<<"#time[fs]\tnu_ee[1/s]\tnu_ei[1/s]\tnu_tot[1/s]\tepsilon_r\tepsilon_i\tintraband_percentage\tT_e[K]\tT_i[K]\n";
        for (unsigned int j = 0; j < n_size; j++) {

 			e_ph_nu[j] =  drude_point->ElectronPhononCollisionFreq(T_lattice[j])/hbar_eV;

			spd_nu[j] = drude_point->get_el_el_const() * n_d[j]  * (n_d[0] - n_d[j]) * std::pow(unit_volume_SI,2.0);
			//spd_nu[j] = drude_point->ElectronElectronCollisionFreq_neq(n_sp[j])/hbar_eV;

			tot_nu[j] = spd_nu[j] + e_ph_nu[j];

			complexDielectricFunction[j] = drude_point->ComplexDielectricDrude_neq(GetOmegaProbe(), T_lattice[j], n_sp[j], n_d[j])
            						+lorentz_point->ComplexDielectricLorentz_neq(GetOmegaProbe(), T_electron[j], mu_sp[j], mu_d[j], n_d[j]);

			double intraband_percentage = fermi_ptr->IntrabandPercentageTotal(T_electron[j]*K_in_au, mu_sp[j]*eV_to_Hartree, mu_d[j]*eV_to_Hartree, GetOmegaProbe()*eV_to_Hartree);

			optics<<Time[j]<<"\t"<<spd_nu[j]<<"\t"<<e_ph_nu[j]<<"\t"<<tot_nu[j]<<"\t"<<complexDielectricFunction[j].real()<<"\t"<<complexDielectricFunction[j].imag()<<"\t"<<intraband_percentage<<"\t"<<T_electron[j]<<"\t"<<T_lattice[j]<<"\n";

			//	n_complex[j] = std::sqrt(complexDielectricFunction[j]);

			//R_s_polarized_film[j] = drude_point->sPolarizedReflectivityFilm(GetOmegaProbe(),complexDielectricFunction[j]);

			//R_p_polarized_film[j] = drude_point->pPolarizedReflectivityFilm(GetOmegaProbe(),complexDielectricFunction[j]);

			//R_n_film[j] = drude_point->NormalIncidenceReflectivityFilm(GetOmegaProbe(),complexDielectricFunction[j]);

			//T_s_polarized_film[j] = drude_point->sPolarizedTransmissivityFilm(GetOmegaProbe(),complexDielectricFunction[j]);

			//T_p_polarized_film[j] = drude_point->pPolarizedTransmissivityFilm(GetOmegaProbe(),complexDielectricFunction[j]);

			//A_p_polarized_film[j] = 1.0 - (R_p_polarized_film[j] + T_p_polarized_film[j]);

			//A_s_polarized_film[j] = 1.0 - (R_s_polarized_film[j] + T_s_polarized_film[j]);

			//sigma_r[j] = vacuum_permittivity * (GetOmegaProbe()/hbar_eV) * complexDielectricFunction[j].imag() * 9e9;  // converted in Hz

			//sigma_i[j] = vacuum_permittivity * (GetOmegaProbe()/hbar_eV) * (drude_point->GetDielectricConst() - complexDielectricFunction[j].real()) * 9e9; //converted in Hz

			//dcConductivity[j] = drude_point->dcConductivity_neq(T_lattice[j], n_sp[j], n_d[j]); // multiply by 9e9 to convert S/m to 1/s


			//write optics you wish to file
			//optics<<Time[j]<<"\t"<<R_p_polarized_film[j]<<"\n";

        }

    }else{
        std::cerr << "\n***The file was not correctly opened. Please check the path where you would like to save the file.***\n";
        exit(1);
    }
    optics.close();
		*/

	/*
	//std::string additional_appendix = "";
	std::string additional_appendix = "_fixedFrequency";
	std::string transientOptics = "var/optics/1dReflec_Trans"+appendix+additional_appendix+"_"+std::to_string((int)(probe_wavelength_const*1e9))+".OUT";
	std::ofstream optics(transientOptics);
	if (optics.is_open()) {
		optics<<std::scientific<<std::setprecision(15);

		optics<<"#omega_probe[eV]\tT\tdelT/T\tR\tdelR/R\n";

		double measurement_time = 7.2; //ps
		unsigned int j = (unsigned int) (measurement_time*1e-12 - simulation_t_start)/simulation_t_step;

		for(unsigned int l = 0; l < probe_wavelength_amount; l++) {
			double omega = min_probe_wavelength + ((max_probe_wavelength-min_probe_wavelength)/(probe_wavelength_amount-1)) * l;

			std::complex<double> complexDielectricFunction_eq = drude_point->ComplexDielectricDrude_neq(omega, 300, 1/unit_volume_SI, 10/unit_volume_SI)
            						+lorentz_point->ComplexDielectricLorentz_neq(omega, 300, fermi_ptr->GetFermiEnergy()*Hartree_to_eV, fermi_ptr->GetFermiEnergy()*Hartree_to_eV, 10/unit_volume_SI);

			std::complex<double> complexDielectricFunction_neq = drude_point->ComplexDielectricDrude_neq(omega, T_lattice[j], n_sp[j], n_d[j])
            				+lorentz_point->ComplexDielectricLorentz_neq(omega, T_electron[j], mu_sp[j], mu_d[j], n_d[j]);

			double R_film = drude_point->pPolarizedReflectivityFilm(omega,complexDielectricFunction_neq);
			double R_film_eq = drude_point->pPolarizedReflectivityFilm(omega,complexDielectricFunction_eq);
			double del_R_film = R_film - R_film_eq;

			double T_film = drude_point->pPolarizedTransmissivityFilm(omega,complexDielectricFunction_neq);
			double T_film_eq = drude_point->pPolarizedTransmissivityFilm(omega,complexDielectricFunction_eq);
			double del_T_film = T_film - T_film_eq;

			optics<<omega<<"\t"<<T_film<<"\t"<<(del_T_film/T_film_eq)<<"\t"<<R_film<<"\t"<<(del_R_film/R_film_eq)<<"\n";
		}


    }else{
        std::cerr << "\n***The file was not correctly opened. Please check the path where you would like to save the file.***\n";
        exit(1);
    }
    optics.close();
	*/

	/*
	//std::string additional_appendix = "";
	std::string additional_appendix = "_fixedFrequency";
	std::string transientTransmissivity = "var/optics/transmissivity_copper"+appendix+additional_appendix+".OUT";
	std::string transientReflectivity = "var/optics/reflectivity_copper"+appendix+additional_appendix+".OUT";
	std::ofstream transmissivityOptics(transientTransmissivity);
	std::ofstream reflectivityOptics(transientReflectivity);
	if (transmissivityOptics.is_open() and reflectivityOptics.is_open()) {
		transmissivityOptics<<std::scientific<<std::setprecision(15);
		reflectivityOptics<<std::scientific<<std::setprecision(15);

		for(unsigned int l = 0; l < probe_wavelength_amount; l++) {
			omega[l] = min_probe_wavelength + ((max_probe_wavelength-min_probe_wavelength)/(probe_wavelength_amount-1)) * l;

			//complexDielectricFunction_eq[l] = drude_point->ComplexDielectricDrude_neq(omega[l], 300, 1/unit_volume_SI, 10/unit_volume_SI)
            //				+lorentz_point->ComplexDielectricLorentz_neq(omega[l], 300, fermi_ptr->GetFermiEnergy()*Hartree_to_eV, fermi_ptr->GetFermiEnergy()*Hartree_to_eV, 10/unit_volume_SI);

			complexDielectricFunction_eq[l] = drude_point->ComplexDielectricDrude_neq(omega[l], T_lattice[0], n_sp[0], n_d[0])
            						+lorentz_point->ComplexDielectricLorentz_neq(omega[l], T_electron[0], mu_sp[0], mu_d[0], n_d[0]);

			R_film_eq[l] = drude_point->pPolarizedReflectivityFilm(omega[l],complexDielectricFunction_eq[l]);

			T_film_eq[l] = drude_point->pPolarizedTransmissivityFilm(omega[l],complexDielectricFunction_eq[l]);
		}

		for (unsigned int j = 0; j < n_size; j++) {

			if(j%time_step_factor == 0){

				for(unsigned int l = 0; l < probe_wavelength_amount; l++) {
					std::complex<double> complexDielectricFunction_neq = drude_point->ComplexDielectricDrude_neq(omega[l], T_lattice[j], n_sp[j], n_d[j])
            						+lorentz_point->ComplexDielectricLorentz_neq(omega[l], T_electron[j], mu_sp[j], mu_d[j], n_d[j]);

					double R_film = drude_point->pPolarizedReflectivityFilm(omega[l],complexDielectricFunction_neq);
					double del_R_film = R_film - R_film_eq[l];

					double T_film = drude_point->pPolarizedTransmissivityFilm(omega[l],complexDielectricFunction_neq);
					double del_T_film = T_film - T_film_eq[l];

					transmissivityOptics<<Time[j]*0.001<<"\t"<<omega[l]<<"\t"<<(del_T_film/T_film_eq[l])<<"\n";
					reflectivityOptics<<Time[j]*0.001<<"\t"<<omega[l]<<"\t"<<(del_R_film/R_film_eq[l])<<"\n";
					//transmissivityOptics<<Time[j]*0.001<<"\t"<<omega<<"\t"<<T_film<<"\n";
					//reflectivityOptics<<Time[j]*0.001<<"\t"<<omega<<"\t"<<R_film<<"\n";
				}
			}
		}
	}else{
        std::cerr << "\n***The file was not correctly opened. Please check the path where you would like to save the file.***\n";
        exit(1);
    }
    transmissivityOptics.close();
	reflectivityOptics.close();
	*/


	std::string additional_appendix = "_DFT";
	//std::string additional_appendix = "_easy";
	double T_ph = 1000;
	std::string transientOptics = "var/optics/scatteringEquilibrium"+std::to_string((int)(T_ph))+"K"+additional_appendix+".OUT";
	std::ofstream optics(transientOptics);
	if (optics.is_open()) {
		optics<<std::scientific<<std::setprecision(15);

		optics<<"#T_e[K]\tnu_ee_Fourment[1/s]\tnu_ee_Kirkwood[1/s]\tnu_ei[1/s]\tnu_c[1/s]\tnu_tot_Fourment[1/s]\tnu_tot_Kirkwood[1/s]\n";

		unsigned int temperature_amount = 10000;
		double min_temperature = 300;
		double max_temperature = 10000;

		double fermi_temperature = 8.16e4;
		double nu_ei = drude_point->ElectronPhononCollisionFreq(T_ph)/hbar_eV;

		for(unsigned int l = 0; l < temperature_amount; l++) {
			double temperature = min_temperature + ((max_temperature-min_temperature)/(temperature_amount-1)) * l;

			double nu_ee_Fourment = drude_point->ElectronElectronCollisionFreq(temperature)/hbar_eV;

			double nu_ee_Kirkwood = (3*k_B/(hbar*fermi_temperature)) * temperature * temperature;

			double nu_c = std::sqrt((2*k_B*fermi_temperature/m_e)+(k_B*temperature/m_e))/91e-12;

			double nu_tot_Fourment = drude_point->TotalCollisionFreq(temperature, T_ph)/hbar_eV;

			double nu_tot_Kirkwood = std::min(nu_ei+nu_ee_Kirkwood, nu_c);

			optics<<temperature<<"\t"<<nu_ee_Fourment<<"\t"<<nu_ee_Kirkwood<<"\t"<<nu_ei<<"\t"<<nu_c<<"\t"<<nu_tot_Fourment<<"\t"<<nu_tot_Kirkwood<<"\n";
		}


    }else{
        std::cerr << "\n***The file was not correctly opened. Please check the path where you would like to save the file.***\n";
        exit(1);
    }
    optics.close();


}
