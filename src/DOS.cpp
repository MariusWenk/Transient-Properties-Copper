/* 
 * File:   used_data.h
 * Author: Wenk, Ag Rethfeld
 *
 * Created on 16. May 2022, 11:21
 */

#include "DOS.h"

DOS::DOS(){
	
}

DOS::DOS(std::string dos_filename){
	integrator_ptr = new Cintegration;

    	if (dos_filename != "") {

        	std::ifstream data(dos_filename);

        	std::vector<double> value(4);
        
		if (data.is_open()) {

            		while (data >> value[0] >> value[1] >> value[2] >> value[3]) {

                		if (value[3] != 0.) {

                    			if (dosfile.size() > 0) {
                        			
						value[0] -= dosfile[0][0]; //fermi energy is 0 in data
                        
						value[0] *= Hartree_to_eV;
                        
						value[1] /= Hartree_to_eV;
                        
						value[2] /= Hartree_to_eV;
                        
						value[3] /= Hartree_to_eV;
                        			dosfile.push_back(value);
                    			}else{
	
        			                dosfile.push_back(value);
                    			}
                		}
            		}
            		fermi_energy = (-1.0 * dosfile[0][0])  * Hartree_to_eV; // directly obtained from input Seb_DOS which had FermiEnergy at 0
            
			std::cout << "\nIn class \"BoostRootFinding1D\", Fermi energy located at " << fermi_energy << " eV in DOS file!\n";
            
			dosfile[0][0] = 0.0;
         
	   		data.close();
            		
			std::cout << "\nMax DOS energy located at " << dosfile[dosfile.size()-1][0] << " eV in total DOS file\n";
		
		}else{
            		std::cerr << "\nThe file: " << dos_filename << " for the realistic DOS is not correctly opened. \nPlease check the filename or it's path!\n";
            		std::exit(1);
        	}
         
	}
    //    for (unsigned int k = 0; k < (unsigned int) dosfile.size(); k++) {
    //
    //        std::cout<<"energy "<<dosfile[k][0]<<"\t"<<"dos_tot "<<dosfile[k][1]<<"\n";
    //         
    //    }
    
	this -> n = dosfile.size()-1;

    	this->energy = new double [(int) dosfile.size()];
    	
	this->tot_DOS = new double [(int) dosfile.size()];
    	
	this->sp_DOS = new double [(int) dosfile.size()];
    	
	this->d_DOS = new double [(int) dosfile.size()];

	this ->distance_d_band_E_F = distance_d_band_E_F_const; // in eV

	this ->sp_band_mass = sp_band_mass_const * m_e;

	this ->d_band_mass = d_band_mass_const * m_e;	
 
	for (unsigned int i = 0; i < (unsigned int) dosfile.size(); i++) {
		energy[i] = dosfile[i][0];
		sp_DOS[i] = dosfile[i][1] * sp_DOS_stretching;
		d_DOS[i] = dosfile[i][2] * d_DOS_stretching; 
		//d_DOS[i] = dosfile[i][2] * SmoothFermiFunction(energy[i]) * d_DOS_stretching;  // d_DOS weighted  to obtain 10d-electrons  
		tot_DOS[i] = sp_DOS[i] + d_DOS[i];
		
	}
	
	total_dos_accelerator = gsl_interp_accel_alloc();
	
	total_dos_spline = gsl_spline_alloc(gsl_interp_akima, (int) dosfile.size());
	gsl_spline_init(total_dos_spline, energy, tot_DOS, (int) dosfile.size());
	
	sp_accel = gsl_interp_accel_alloc();
	sp_spline = gsl_spline_alloc(gsl_interp_akima, (int) dosfile.size());
	
	gsl_spline_init(sp_spline, energy, sp_DOS, (int) dosfile.size());
	d_accel = gsl_interp_accel_alloc();
	
	d_spline = gsl_spline_alloc(gsl_interp_akima, (int) dosfile.size());
	gsl_spline_init(d_spline, energy, d_DOS, (int) dosfile.size());

	this ->temperature_start = temperature_start_const;

	this ->LowerBound = energy[0];

	this ->UpperBound = energy[dosfile.size()-2]; 
   
	//this ->UpperBound = dosfile[dosfile.size() - 1] [0];

	this ->unitVolume = unit_volume_au / std::pow(m_to_a0, 3);
}

DOS::~DOS(){

}

void read_DOS(std::string dos_filename){
	
    
}