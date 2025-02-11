

/******************************
 * File:   main.cpp
 * Author: Ndione @ Ag Rethfeld
 ******************************
 ******************************/

#include "main.h"

int main(int argc, char** argv) {

	/*MPI_Init(NULL, NULL);

	int world_size;
	int world_rank;

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	*/

	std::cout.precision(15);
	std::cout.setf(std::ios::scientific);

	std::clock_t start_time = std::clock();
	time_t t = time(NULL);
 	tm* timePtr = localtime(&t);

	int index;

	std::cout<<"\nTo write data in a file, please choose the right integer\n"
		<<"0 saves data of 2band1T\n"
		<<"1 saves data of 2band2T\n"
		<<"2 saves data of 3band1T\n"
		<<"3 saves data of 3band2T\n"
		<<"4 saves data of 4band1T\n"
		<<"5 saves data of 4band2T\n"
		<<"Any other key to exit\n";

	std::cin >> index;

	std::cout<<"\nDynamics started. ODE solver running\n";

	while (std::cin.fail()){

		std::cerr<<"\n\33[1;31mFailed. Unknow type. Use an interger.\33[0m\n";
		std::cin.clear();
		std::exit(1);
	}

	std::cout<<"\nYou chose index idx = "<<index<<"\n";

	solver solver(index);

	std::cout<<"\nDynamics ended\n";


	std::cout<<"*****************************************"
		<<"\nTransient optical properties starts here\n";

	BoostRootFinding1D boost_obj(DOS_input);
	dielectric_function base_obj;
	dielectric_function_drude drude_obj(&boost_obj);
	dielectric_function_lorentz lorentz_obj(&boost_obj);

	std::string appendix = file_name_appendix;
	std::string inputFile = "var/FullDynamics"+appendix+".in";

 	transient_optical_param param_obj(inputFile, &drude_obj, &lorentz_obj);
	//param_obj.OpticalPropertiesDrudeLorentz();

	std::cout<<"\nTransient optical properties ended here\n";

	std::cout<<"*****************************************"
		<<"\nAbsorption properties starts here\n";

	//Cfermi_distribution *fermi_ptr = new Cfermi_distribution("input/DOS_Cu/Cu_l_proj_dos.in");
	Cfermi_distribution *fermi_ptr = new Cfermi_distribution("input/DOS_Cu/DOS_easy_Cu.in");
	std::string absorptionFile = "var/absorption_ratio.OUT";
	std::ofstream absorption(absorptionFile);
	if (absorption.is_open()) {
		absorption<<std::scientific<<std::setprecision(15);
		double max_photo_energy = 60;
		double min_photo_energy = 0;
		double T1 = 2000;
		unsigned int photo_energy_steps = 1000;
		for(unsigned int j = 0; j < photo_energy_steps; j++) {
			double energy = min_photo_energy + ((max_photo_energy-min_photo_energy)/(photo_energy_steps-1)) * j;
			double interband_percentage_T0 = fermi_ptr->InterbandPercentageTotal(300*K_in_au, fermi_ptr->GetFermiEnergy(), fermi_ptr->GetFermiEnergy(), energy * eV_to_Hartree); //All energies in Cfermi_distribution in Hartrees
			double intraband_percentage_T0 = fermi_ptr->IntrabandPercentageTotal(300*K_in_au, fermi_ptr->GetFermiEnergy(), fermi_ptr->GetFermiEnergy(), energy * eV_to_Hartree);
			double interband_percentage_T1 = fermi_ptr->InterbandPercentageTotal(T1*K_in_au, fermi_ptr->GetFermiEnergy(), fermi_ptr->GetFermiEnergy(), energy * eV_to_Hartree);
			double intraband_percentage_T1 = fermi_ptr->IntrabandPercentageTotal(T1*K_in_au, fermi_ptr->GetFermiEnergy(), fermi_ptr->GetFermiEnergy(), energy * eV_to_Hartree);

			absorption<<energy<<"\t"<<interband_percentage_T0<<"\t"<<intraband_percentage_T0<<"\t"<<interband_percentage_T1<<"\t"<<intraband_percentage_T1<<"\n";
		}
	}
	absorption.close();

	std::cout<<"\nAbsorption properties ended here\n";

	/*

	std::cout<<"\n*****************************************"
		<<"\nConvolution started here\n";

	std::string outputFile = "var/TransientOptics"+appendix+"_"+std::to_string((int)(probe_wavelength_const*1e9))+".OUT";

	convolution conv_obj(outputFile);

	std::string convolvedOptics = "var/transientOpticsConvolved.OUT";
	std::ofstream file(convolvedOptics);
	if(file.is_open()){
		double time = -1e-12;
		while(time<10e-12){

			file<<time*1e15
			<<"\t"<<conv_obj.ConvolutedReflectivity_pPol(time, -2e-12 , 7e-12)
			<<"\n";
			time+=1e-16;
		}
	}else{
		std::cerr<<"\nFile "<<convolvedOptics<<" is not opened!\nExit!\n";
		std::exit(1);
	}
	file.close();

	std::cout<<"\nConvolution ended here\n";

	*/


	start_time = clock()-start_time;


	std::cout
			<<"\n****************************************************************************************************************************"
			<< "\n*DF_Program compiled on "
			<< (timePtr->tm_mday)
			<< ":"
			<< (timePtr->tm_mon)+1
			<< ":"<<(timePtr->tm_year)+1900
		    	<< " "<<"at local time "
			<< (timePtr->tm_hour)
			<< ":"
			<< (timePtr->tm_min)
			<< ":"
			<< (timePtr->tm_sec)
			<< "\n"

			<< "\n*The real time: <<Time between the instant that the process was spawned until it ended>> is: "
			<< (double)start_time/CLOCKS_PER_SEC/60.0
			<<" minuts\n"

			<< "\n*The user time: <<Time spent on the process by the cpu>>: is: ....do it later "
			<< "\n"

			<< "\n*The sytem time: <<Time spent in the OS while performing actions requested by the process>> is:...do it later "
			<<"\n"
																																       				<<"\n*****************************************************************************************************************************\n";

	// MPI_Finalize();

	std::cout<<"\nDF_Program over.\n";

	return EXIT_SUCCESS;

}
