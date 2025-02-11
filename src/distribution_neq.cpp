/************************************ 
 * File:   fermi_distribution.h
 * Author: Ndione, Ag Rethfeld
 ***********************************
 * Created on 09. January 2019, 10:55
 ************************************/

#include "distribution_neq.h"


Cdistribution_neq::Cdistribution_neq(){
	std::cerr<<"\nUse other constructor\n exits\n";
	std::exit(1);
}
   
Cdistribution_neq::Cdistribution_neq(BoostRootFinding1D *boost_pointer, std::string dos_file, double T_elec, double mu_sp, double mu_d){
	
		this->T_elec = T_elec;
		this->mu_sp = mu_sp;
		this->mu_d = mu_d;

	 	if (dos_file != "") {

        	std::ifstream data(dos_file);

        	std::vector<double> value(4);

		if (data.is_open()) {

            		while (data >> value[0] >> value[1] >> value[2] >> value[3]) {

                		if (value[3] != 0.) {

                    			if (load.size() > 0) {    
           
						value[0] -= load[0][0];   
                  
						value[0] *= Hartree_to_eV;

						value[1] /= Hartree_to_eV;

						value[2] /= Hartree_to_eV; 
       
						value[3] /= Hartree_to_eV;

                        			load.push_back(value);
                    			}else{
        			                load.push_back(value);
                    			}
                		}
            		}            	
            
			load[0][0] = 0.0;

	   		data.close();
            		
		}else{
            		std::cerr << "\nThe file: " << dos_file << " for the realistic DOS is not correctly opened. \nPlease check the filename or it's path!\n";
            		std::exit(1);
        	}
         
	}

	this->E = new double [(int) load.size()]; 
  	
	this->totDOS = new double [(int) load.size()];  
	
	this->spDOS = new double [(int) load.size()]; 
	
	this->dDOS = new double [(int) load.size()];

	for (unsigned int i = 0; i < (unsigned int) load.size(); i++) {		
		E[i] = load[i][0];

		spDOS[i] = load[i][1] * sp_DOS_stretching;

		dDOS[i] = load[i][2] * d_DOS_stretching; 
		 
		totDOS[i] = spDOS[i] + dDOS[i];


	}	

	this ->boost_pointer = boost_pointer;

	integration_pointer = new Cintegration;
}     
               
Cdistribution_neq::~Cdistribution_neq(){

	delete integration_pointer;	 
}

double Cdistribution_neq::GetTemperature(){

	return T_elec;
}

double Cdistribution_neq::GetMu_sp(){

	return mu_sp;
}

double Cdistribution_neq::GetMu_d(){

	return mu_d;
}
 
double Cdistribution_neq::GetDistribution_sp(double energy){

	return boost_pointer->GetFermiDistribution(energy, GetTemperature(), GetMu_sp());
}

double Cdistribution_neq::GetDistribution_d(double energy){

	return boost_pointer->GetFermiDistribution(energy, GetTemperature(), GetMu_d());
}

double Cdistribution_neq::GetPhiFunction_sp(double energy){

	return -std::log((1.0/GetDistribution_sp(energy)) - 1.0);
}

double Cdistribution_neq::GetPhiFunction_d(double energy){

	return -std::log((1.0/GetDistribution_d(energy)) - 1.0);
}

double Cdistribution_neq::NumberDensity_sp(double energy){
	
	return GetDistribution_sp(energy) * boost_pointer->GetDOS_sp(energy);
}

double Cdistribution_neq::NumberDensity_d(double energy){

	return GetDistribution_d(energy) * boost_pointer->GetDOS_d(energy);
}

double Cdistribution_neq::GetDistributionNonEq(double energy){

	return (NumberDensity_sp(energy) + NumberDensity_d(energy)) / (boost_pointer->GetDOS_sp(energy) +  boost_pointer->GetDOS_d(energy));
}

double Cdistribution_neq::GetPhiFunctionNonEq(double energy){
	
	return -std::log((1.0/GetDistributionNonEq(energy)) - 1.0);
} 

double Cdistribution_neq::ParticlesFromDistributionNonEq(){

	auto integrand  =  [ = ] (double energy)->double {
		
		return GetDistributionNonEq(energy) * boost_pointer->GetDOS(energy);
	};

	return integration_pointer->integrate_qag(integrand, boost_pointer->GetLowerBound(), boost_pointer->GetUpperBound());
	  	 
}

double Cdistribution_neq::ParticlesFromFermiDistribution(){

	 return boost_pointer->NumberofParticles(boost_pointer->GetTemperatureStart(), boost_pointer->GetFermiEnergy());

}

void Cdistribution_neq::eq_distribution(std::string filename){

	std::ofstream output_file(filename);
	if(output_file.is_open()){
		output_file<< std::scientific << std::setprecision(15);
		double root_mu = boost_pointer->GetMuRootFinding(GetTemperature());
		output_file<<"#Te = "<<GetTemperature()<<" [K]"
		<<"\n#mu_eq-E_F = "<<root_mu-boost_pointer->GetFermiEnergy()<<" [eV]"
		<<"\n#energy-E_F [eV]\teq_Fermi_distribution\n";
		for (unsigned int i=0;  i< load.size()-1; i++){
			output_file<<E[i]-boost_pointer->GetFermiEnergy()
			<<"\t"<<boost_pointer->GetFermiDistribution(E[i], GetTemperature(), root_mu)
			<<"\n";
		}
	}else{
		std::cerr<<"\nFile "<<filename<<" is not opened!\n";
		std::exit(1);
	}
	output_file.close();
}

void Cdistribution_neq::partialDistribution(std::string filename){

	std::ofstream output_file(filename);
	if(output_file.is_open()){
		output_file<< std::scientific << std::setprecision(15);
		output_file<<"#Te = "<<GetTemperature()<<" [K]"
		<<"\n#mu_sp-E_F = "<<GetMu_sp()-boost_pointer->GetFermiEnergy()<<" [eV]\n"
		<<"\n#mu_d-E_F = "<<GetMu_d()-boost_pointer->GetFermiEnergy()<<" [eV]\n"
		<<"\n#energy-E_F [eV]\tf_s\tf_d\n";
		for (unsigned int i=0;  i< load.size()-1; i++){
			output_file<<E[i]-boost_pointer->GetFermiEnergy()
				<<"\t"<<GetDistribution_sp(E[i])
				<<"\t"<<GetDistribution_d(E[i])
				<<"\n";
		}
	}else{
		std::cerr<<"\nFile "<<filename<<" is not opened!\n";
		std::exit(1);
	}
	output_file.close();
}
void Cdistribution_neq::GetDistributionNonEq(std::string filename){

	std::ofstream output_file(filename);
	double *f_sp = new double [load.size()];
	double *f_d  = new double [load.size()];
	if(output_file.is_open()){
		output_file<< std::scientific << std::setprecision(15);
		output_file<<"#quasi_Te = "<<GetTemperature()<<" [K]"
		<<"\n#mu_sp-E_F = "<<GetMu_sp()-boost_pointer->GetFermiEnergy()<<" eV"
		<<"\n#mu_d-E_F = "<<GetMu_d()-boost_pointer->GetFermiEnergy()<<" eV"
		<<"\n#energy-E_F [eV]\tNEDF\n";
		for (unsigned int k=0; k < load.size()-1 ; k++){
			f_sp[k] = GetDistribution_sp(E[k]);
			f_d[k]  = GetDistribution_d(E[k]);
		
			output_file<<E[k]-boost_pointer->GetFermiEnergy()
			<<"\t"
			<<((f_sp[k]*spDOS[k] + f_d[k]*dDOS[k]) / (spDOS[k] + dDOS[k]))	
			<<"\n";
		}
	}else{
		std::cerr<<"\nFile "<<filename<<" is not opened!\n";
		std::exit(1);
	}
	output_file.close();	
	delete f_sp;
	delete f_d;
}














