#include "nonThermalOptics.hpp"



nonThermalOptics::nonThermalOptics(){
	std::cerr<<"\nUse the other constructors\n";
	std::exit(1);
}

nonThermalOptics::nonThermalOptics(std::string non_thermal_distribution){

	if(boost::filesystem::exists(non_thermal_distribution)){
	
		std::ifstream readData(non_thermal_distribution);
		
		std::vector<double> numColumn(2);
		
		if(readData.is_open()){ // replace this line later
		
			while (readData >> numColumn[0] >> numColumn[1]){
			
				numColumn[0] *= 1; //Hartree_to_eV; no need to convert energy 
				numColumn[1] *= 1; // unitless distribution
				fileLoader.push_back(numColumn); 	
			}	
			
		}else{
			std::cerr<<"\nFile exists but is empty\n";
			std::exit(1);
		}
		
		readData.close();
		
		this->energy = new double [fileLoader.size()];
		
		this->distribution_f = new double [fileLoader.size()];	

		for (int id = 0; id != (int) fileLoader.size(); id++){
		
			energy[id]         = fileLoader[id][0];
			distribution_f[id] = fileLoader[id][1];
			//std::cout<<"index "<<id<<"\tenergy "<<energy[id]<<"\tdistribution "<<distribution_f[id]<<"\n";
		}
		
		
		f_accel = gsl_interp_accel_alloc();

		f_spline = gsl_spline_alloc(gsl_interp_akima, fileLoader.size());

		gsl_spline_init(f_spline, energy, distribution_f, fileLoader.size());

	}else{
		std::cerr<<"\n\33[1;31mFilename "<<non_thermal_distribution
		<<" not found. Check the path or the file's name. \nThe program exits. \33[0m\n";
		std::exit(1);
	}
	
}

nonThermalOptics::~nonThermalOptics(){

	delete [] energy;
	delete [] distribution_f;
	
	gsl_spline_free(f_spline);
	gsl_interp_accel_free(f_accel);

}


double nonThermalOptics::interpolateDistribution(double Energy){

	if((Energy < energy[0]) || Energy > energy[(int) fileLoader.size()-2]){
		std::cerr<<"\n\33[1;31mEnery value out of range \33[0m\n";
		std::exit(1);
	}
	
	return gsl_spline_eval(f_spline, Energy, f_accel);
}
























