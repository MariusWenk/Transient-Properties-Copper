/* 
 * File:   fermi_distribution.cpp
 * Author: ndione, Ag Rethfeld
 * 
 * Created on 22. Februar 2018, 09:48  
 */

#include "fermi_distribution.h"

Cfermi_distribution::Cfermi_distribution(std::string dos_filename) {
	
    	this ->dos_factor = 2; // Should not multiply the DOS by dosfactor since already done

	sp_band_pointer = new Sp_band;

    	d_band_pointer = new D_band;

    	integrator_ptr = new Cintegration;

	this ->temperature_start = temperature_start_const * K_in_au;

	if (dos_filename != "") {

        	std::ifstream data(dos_filename);

        	std::vector<double> value(4);
        
		if (data.is_open()) {

            		while (data >> value[0] >> value[1] >> value[2] >> value[3]) {

                		if (value[3] != 0.) {//only take into account data where the total dos is not zero

                    			//convert states/Hartree/atom to 1/Hartree/a_0³
                    			value[1] /= (unit_volume_au);

                    			value[2] /= (unit_volume_au);

                    			value[3] /= (unit_volume_au);

                    				if (dosfile.size() > 0) {
                        				
							value[0] -= dosfile[0][0]; //fermi energy is 0 in data. Like this, fermi energy is -dosfile[0][0]
                        				value[0] *= 1.0; // The energy is already in Hartree, we don't need a conversion
                       					//value[1] *= 1.0; 
                        				//value[2] *= 1.0;
                        				//value[3] *= 1.0;

                        				dosfile.push_back(value);

                    				}else{
                        			
							dosfile.push_back(value);						
						}
				}

	
			}
           
			fermi_energy = -1.0 * dosfile[0][0];
	    		//fermi_energy = (-1.0 * dosfile[0][0]) + 0.29266027; // fermi energy previously located at 0.29266027 Hartree in Nils DOS
            
			std::cout << "\nIn class \"fermi_distribution\",Fermi energy located at " << fermi_energy << " Hartrees = " << fermi_energy * Hartree_to_eV << " eV in DOS file!\n";

            		dosfile[0][0] = 0.0;

            		data.close();

            		std::cout << "\nMax DOS energy located at " << dosfile[dosfile.size() - 1][0] << " Hartrees = " << dosfile[dosfile.size() - 1][0] * Hartree_to_eV << " eV\n";

		}else{
            		std::cerr << "\nThe file: " << dos_filename << " for the realistic DOS is not correctly opened. \nPlease check the filename or it's path!\n";

            		std::exit(1);
        	}

		data.close();

	}
	
	this ->unit_volume = unit_volume_au;

    	this ->n = (unsigned int) dosfile.size()-1; // number of discretization points of the DOS file
    	
	std::cout<<"\nThe number of discretization points of the DOS file is  n = "<<n<<"\n";
    	
	this->energy = new double [(int) dosfile.size()];
    	
	this->DOS_sp = new double [(int) dosfile.size()];
    	
	this->DOS_d = new double [(int) dosfile.size()];
    	
	this->DOS = new double [(int) dosfile.size()];

	this ->distance_d_band_E_F = distance_d_band_E_F_const * eV_to_Hartree;
 
 	for (unsigned int i = 0; i != dosfile.size(); i++) {
		
        	energy[i] = dosfile[i][0];

		DOS_sp[i] = dosfile[i][1] * sp_DOS_stretching;
 
		DOS_d[i] = dosfile[i][2] * d_DOS_stretching;
		//DOS_d[i] = dosfile[i][2] * SmoothFermiFunction(energy[i]) * d_DOS_stretching;
		 
		DOS[i] = DOS_sp[i] + DOS_d[i];
    	}
   
	sp_dos_accel = gsl_interp_accel_alloc();

	sp_dos_spline = gsl_spline_alloc(gsl_interp_akima, (int) dosfile.size());

	gsl_spline_init(sp_dos_spline, energy, DOS_sp, (int) dosfile.size());

	d_dos_accel = gsl_interp_accel_alloc();

	d_dos_spline = gsl_spline_alloc(gsl_interp_akima, (int) dosfile.size());

	gsl_spline_init(d_dos_spline, energy, DOS_d, (int) dosfile.size());

	tot_dos_accel = gsl_interp_accel_alloc();

	tot_dos_spline = gsl_spline_alloc(gsl_interp_akima, (int) dosfile.size());

	gsl_spline_init(tot_dos_spline, energy, DOS, (int) dosfile.size());
}

Cfermi_distribution::Cfermi_distribution() {

	sp_band_pointer = new Sp_band;

    	d_band_pointer = new D_band;

	integrator_ptr = new Cintegration;

    	this ->fermi_energy = fermi_energy_const * eV_to_Hartree;
    
	this ->dos_factor = 2;

    	this ->n = 4000;

    	energy = new double [n];

    	DOS = new double [n];

    	DOS_sp = new double [n];

    	DOS_d = new double [n];

    	double delta_E = 2.0 * fermi_energy / n;

    	for (unsigned int i = 1; i < n; i++) {

        	energy[i] = i*delta_E;
        
		DOS[i] = (sp_band_pointer->Get_sp_DOS(energy[i]) + d_band_pointer->Get_d_DOS(energy[i]))  ;

	        DOS_sp[i] = (sp_band_pointer->Get_sp_DOS(energy[i]))  ;
        
		DOS_d[i] = d_band_pointer->Get_d_DOS(energy[i])  ;
		//std::cout<<energy[i]*Hartree_to_eV<<"\t"<<DOS_sp[i]/(Hartree_to_Joule*std::pow(a0_to_m, 3.))<<"\t"<<DOS_d[i]<<"\t"<<DOS[i]<<"\n";
    }
}
Cfermi_distribution::~Cfermi_distribution() {

	delete [] energy;

	delete [] DOS;

 	delete [] DOS_sp;

	delete [] DOS_d;

        delete integrator_ptr;

	delete sp_band_pointer;

	delete d_band_pointer;

	gsl_spline_free(sp_dos_spline);

	gsl_interp_accel_free(sp_dos_accel);

	gsl_spline_free(d_dos_spline);

	gsl_interp_accel_free(d_dos_accel);

	gsl_spline_free(tot_dos_spline);

	gsl_interp_accel_free(tot_dos_accel);

}

double Cfermi_distribution::GetUnitVolume(){

	return unit_volume;
}

double Cfermi_distribution::GetDistance_d_BandtoFermi(){

	return distance_d_band_E_F;
}

double Cfermi_distribution::GetFermiEnergy() {

    	return fermi_energy;
}
double Cfermi_distribution::GetTemperatureStart(){

	return temperature_start;
}
double Cfermi_distribution::LowerBound() {

    	return energy[0];
    	//return dosfile[0][0];
}

double Cfermi_distribution::UpperBound() {

    	return energy[dosfile.size() - 1];
    //    return dosfile[dosfile.size()-1][0];
}

double Cfermi_distribution::GetNumberDiscPoints() {

    	return n;
}

double Cfermi_distribution::StepFunction(double energy){
	
	if ((energy >= 0.0) && (energy < (GetFermiEnergy() -  GetDistance_d_BandtoFermi()))){

		return  1.0;
	}else{  
	
		return 0.0;		
	}
}

double Cfermi_distribution::GetFermiDistribution(double energy, double T, double mu) {

	return 1.0 / (exp((energy - mu) / (k_B_au * T)) + 1.0);

}

double Cfermi_distribution::SmoothFermiFunction(double energy){

	double temp = temperature_smooth_fermi_function * K_in_au;
	
	//double mu = GetFermiEnergy() - GetDistance_d_BandtoFermi();
	double mu = GetFermiEnergy() - 1*eV_to_Hartree;

	return GetFermiDistribution(energy, temp, mu);
	
}

double Cfermi_distribution::GetDOS_sp(double energy) {

	if (dosfile.size() > 0) {

        	if ((energy < dosfile[0][0]) || (energy > dosfile[dosfile.size() - 2][0])) {

            		return 0.;

        	}else{
		 
			return gsl_spline_eval(sp_dos_spline, energy, sp_dos_accel);
        
    		} 
	}else {

       		 std::cerr << "\nNo points found in sp-DOS that you want to interpolate! The programm exits.\n";
       		 std::exit(1);
  	 }
}

double Cfermi_distribution::GetDOS_d(double energy) {

	 if (dosfile.size() > 0) {

        	if ((energy < dosfile[0][0]) || (energy > dosfile[dosfile.size() - 2][0])) {

            		return 0.;
       		 }else{
					 
			return gsl_spline_eval(d_dos_spline, energy, d_dos_accel) * SmoothFermiFunction(energy);
     		}
   	 } else {

        	 std::cerr << "\nNo points found in d-DOS that you want to interpolate! The programm exits.\n";
       		 std::exit(1);
   	 }
}

double Cfermi_distribution::GetDOSTot(double energy) {
    
	if (dosfile.size() > 0) {

        	if ((energy < dosfile[0][0]) || (energy > dosfile[dosfile.size() - 2][0])) {

            		return 0.;
       		 }else{
		 
			return gsl_spline_eval(tot_dos_spline, energy, tot_dos_accel) * SmoothFermiFunction(energy);
		}
	}else{

       		 std::cerr << "\nNo points found in tot-DOS that you want to interpolate! The programm exits.\n";
       		 std::exit(1);
       }
}

double Cfermi_distribution::NumberStatesTotal_sp(double T, double mu_sp, double mu_d, double photon_energy) {
	 /* if (T > .0) {
        
         auto integrand = [ = ](double Energy) -> double {

             return GetDOS_sp(Energy) * GetFermiDistribution(Energy, T, mu)*(1.0 - GetFermiDistribution(Energy + photon_energy, T, mu));

         };

         return integrator_ptr->integrate_qag(integrand, LowerBound(), UpperBound());
     } else {

         std::cerr << "\nPlease start calculation at temperature greater than zero!\n";
        
         std::exit(1);
     }*/
    	double result = 0.0;

		if (T>0.){

    			for (unsigned int j = 1; j < n-1; j++) {

        			result += DOS_sp[j] * GetFermiDistribution(energy[j], T, mu_sp) * (1.0 - GetFermiDistribution(energy[j] + photon_energy, T, mu_sp)) * (energy[j + 1] - energy[j - 1]);

    			}

    			return 0.5 * (result

            		+ DOS_sp[0] * GetFermiDistribution(energy[0], T, mu_sp) * (1.0 - GetFermiDistribution(energy[0] + photon_energy, T, mu_sp)) * (energy[1] - energy[0])
            		+ DOS_sp[n - 1] * GetFermiDistribution(energy[n - 1], T, mu_sp) * (1.0 - GetFermiDistribution(energy[n - 1] + photon_energy, T, mu_sp)) * (energy[n - 1] - energy[n - 2])

            		);
		}else{

		for (unsigned int k=1; k < n-1; k++){
	
			result += DOS_sp[k]*GetFermiDistribution(energy[k], T, mu_sp)*(1.0-GetFermiDistribution(energy[k]+photon_energy, T, mu_sp))*(energy[k+1] - energy[k-1]);
		}
	
		return 0.5 * (result
            			 +DOS_sp[0] * GetFermiDistribution(energy[0], T, mu_sp) * (1.0 - GetFermiDistribution(energy[0] + photon_energy, T, mu_sp)) * (energy[1] - energy[0])
           			 +DOS_sp[n - 1] * GetFermiDistribution(energy[n - 1], T, mu_sp) * (1.0 - GetFermiDistribution(energy[n - 1] + photon_energy, T, mu_sp)) * (energy[n - 1] - energy[n - 2])
            			);
		
		}
}

double Cfermi_distribution::NumberStatesTotal_d(double T, double mu_sp, double mu_d, double photon_energy) {

	/*if (T > .0) {
        
         auto integrand = [ = ](double Energy) -> double {

             return GetDOS_d(Energy) * GetFermiDistribution(Energy, T, mu)*(1.0 - GetFermiDistribution(Energy + photon_energy, T, mu));

         };

         return integrator_ptr->integrate_qag(integrand, LowerBound(), UpperBound());
     } else {

         std::cerr << "\nPlease start calculation at temperature greater than zero!\n";

         std::exit(1);
     }*/
    	double result = 0.0;
	
	if (T>0.){
   
	 	for (unsigned int j = 1; j < n - 1; j++) {

        		result += DOS_d[j] * GetFermiDistribution(energy[j], T, mu_d) * (1.0 - GetFermiDistribution(energy[j] + photon_energy, T, mu_sp)) * (energy[j + 1] - energy[j - 1]);

         	}

    		return 0.5 * (result

            		+ DOS_d[0] * GetFermiDistribution(energy[0], T, mu_d) * (1.0 - GetFermiDistribution(energy[0] + photon_energy, T, mu_sp)) * (energy[1] - energy[0])
            		+ DOS_d[n - 1] * GetFermiDistribution(energy[n - 1], T, mu_d) * (1.0 - GetFermiDistribution(energy[n - 1] + photon_energy, T, mu_sp)) * (energy[n - 1] - energy[n - 2])

            		);
	
	}else{

		for (unsigned int j=1; j < n-1; j++){
		
			result += (energy[j] > GetFermiEnergy()) ? 0. : DOS_d[j]*GetFermiDistribution(energy[j], T, mu_d)*(1.0-GetFermiDistribution(energy[j]+photon_energy, T, mu_sp))*(energy[j+1]-energy[j-1]);
		}
	
		return 0.5 *(result
				+DOS_d[0]*GetFermiDistribution(energy[0], T, mu_d)*(1.0-GetFermiDistribution(energy[0]+photon_energy, T, mu_sp))*(energy[1]-energy[0])
				+DOS_d[n-1]*GetFermiDistribution(energy[n-1], T, mu_d)*(1.0-GetFermiDistribution(energy[n-1]+photon_energy, T, mu_sp))*(energy[n-1]-energy[n-2])			
			);
	}

}

double Cfermi_distribution::NumberStatesTotal(double T, double mu_sp, double mu_d, double photon_energy){

	return NumberStatesTotal_sp(T, mu_sp, mu_d, photon_energy) + NumberStatesTotal_d(T, mu_sp, mu_d, photon_energy);
}

double Cfermi_distribution::InterbandPercentageTotal(double T, double mu_sp, double mu_d, double photon_energy){

	return NumberStatesTotal_d(T, mu_sp, mu_d, photon_energy)/NumberStatesTotal(T, mu_sp, mu_d, photon_energy);

}

double Cfermi_distribution::IntrabandPercentageTotal(double T, double mu_sp, double mu_d, double photon_energy){

	return NumberStatesTotal_sp(T, mu_sp, mu_d, photon_energy)/NumberStatesTotal(T, mu_sp, mu_d, photon_energy);
}


double Cfermi_distribution::CalcInternalEnergy(double T, double mu) {

    /*   if (T > 0.0) {
        
           auto integrand = [ = ] (double Energy) ->double {

               return Energy * GetDOSTot(Energy) * GetFermiDistribution(Energy, T, mu);
           };
           return integrator_ptr->integrate_qag(integrand, LowerBound(), UpperBound());

       } else {

           std::cerr << "\nPleae start calculation at temperature greater than zero!\n";

           std::exit(1);
       }*/
    	double sum = 0.0;
	
	if (T>0.){
    
		for (unsigned int i = 1; i < n - 1; i++) {

        		sum += (energy[i] * DOS[i] * GetFermiDistribution(energy[i], T, mu))*(energy[i + 1] - energy[i - 1]);
    		}

    		return 0.5 * (sum
            		+ energy[0] * DOS[0] * GetFermiDistribution(energy[0], T, mu)*(energy[1] - energy[0])
            		+ energy[n - 1] * DOS[n - 1] * GetFermiDistribution(energy[n - 1], T, mu)*(energy[n - 1] - energy[n - 2])
            		);
	}else{
		for(unsigned int i=1; i< n-1; i++){
		
			sum += (energy[i] > GetFermiEnergy()) ? 0.0 : DOS[i]*(energy[i+1] - energy[i-1]);
		}

		return 0.5 * (sum
	  	 		+DOS[0]*(energy[1]-energy[0])
				+DOS[n-1]*(energy[n-1] - energy[n-2])
			);
	}
}

double Cfermi_distribution::CalcInternalEnergy_sp(double T, double mu) {

    /*   if (T > .0) {
        
           auto integrand = [ = ] (double Energy) ->double {

               return Energy * GetDOS_sp(Energy) * GetFermiDistribution(Energy, T, mu);
           };

           return integrator_ptr->integrate_qag(integrand, LowerBound(), UpperBound());

       }else {

           std::cerr << "\nPlease start calculation at temperature greater than zero!\n";

           std::exit(1);
       }*/
    	double sum = 0.0;

	if (T>0.0){

    		for (unsigned int i = 1; i < n - 1; i++) {

       			sum += (energy[i] * DOS_sp[i] * GetFermiDistribution(energy[i], T, mu))*(energy[i + 1] - energy[i - 1]);
    		}

    		return 0.5 * (sum
            		+ energy[0] * DOS_sp[0] * GetFermiDistribution(energy[0], T, mu)*(energy[1] - energy[0])
            		+ energy[n - 1] * DOS_sp[n - 1] * GetFermiDistribution(energy[n - 1], T, mu)*(energy[n - 1] - energy[n - 2])
            		);

	}else{
	
		for (unsigned int i = 1; i < n-1; i++){

			sum+= (energy[i] > fermi_energy) ? 0.0 : DOS_sp[i]*(energy[i+1] - energy[i-1]);
		}

		return 0.5 * (
				sum 
				+ DOS_sp[0]*(energy[1] - energy[0])
				+ DOS_sp[n-1]*(energy[n-1] - energy[n-2])
			);
	}
}

double Cfermi_distribution::CalcInternalEnergy_d(double T, double mu) {

    /* if (T > .0) {
        
         auto integrand = [ = ] (double Energy) ->double {

             return Energy * GetDOS_d(Energy) * GetFermiDistribution(Energy, T, mu);
         };

         return integrator_ptr->integrate_qag(integrand, LowerBound(), UpperBound());

     } else {

         std::cerr << "\nPlease start calculation at temperature greater than zero!\n";

         std::exit(1);
     }*/
    	double sum = 0.0;

	if (T>0.){

    		for (unsigned int i = 1; i < n - 1; i++) {

       			sum += (energy[i] * DOS_d[i] * GetFermiDistribution(energy[i], T, mu))*(energy[i + 1] - energy[i - 1]);
    		}

    		return 0.5 * (sum
            		+ energy[0] * DOS_d[0] * GetFermiDistribution(energy[0], T, mu)*(energy[1] - energy[0])
            		+ energy[n - 1] * DOS_d[n - 1] * GetFermiDistribution(energy[n - 1], T, mu)*(energy[n - 1] - energy[n - 2])
			);

	}else{
		for (unsigned int k=1; k< n-1; k++){
			
			sum += (energy[k] > GetFermiEnergy()) ? 0.0 : DOS_d[k]*(energy[k+1]-energy[k-1]);
		}
			
		return 0.5*(sum 
			+ DOS_d[0]*(energy[1]-energy[0])
			+DOS_d[n-1]*(energy[n-1]-energy[n-2])
			);
	}
}

double Cfermi_distribution::CalcDensity(double T, double mu) {

    /* if (T > 0.) {

         auto integrand = [ = ] (double Energy) ->double {

             return GetDOSTot(Energy) * GetFermiDistribution(Energy, T, mu);
         };


         return integrator_ptr->integrate_qag(integrand, LowerBound(), UpperBound());

     } else {

         std::cerr << "\nPlease start calculation at temperature greater than zero!\n";

         std::exit(1);
         }*/
    	double sum = 0.0;

    	if (T > 0.0) {

        	for (unsigned int i = 1; i < n - 1; i++) {

            		sum += (DOS[i] * GetFermiDistribution(energy[i], T, mu))*(energy[i + 1] - energy[i - 1]);
        	}

        	return 0.5 * (sum
                	+ DOS[0] * GetFermiDistribution(energy[0], T, mu)*(energy[1] - energy[0])
                	+ DOS[n - 1] * GetFermiDistribution(energy[n - 1], T, mu)*(energy[n - 1] - energy[n - 2])
                	);
    	} else {

        	for (unsigned int i = 1; i < n - 1; i++) {

            		sum += (energy[i] > fermi_energy) ? 0.0 : DOS[i]*(energy[i + 1] - energy[i - 1]);
        	}

		return 0.5 * (sum 
				+ DOS[0]*(energy[1] - energy[0]) 
				+ DOS[n - 1]*(energy[n - 1] - energy[n - 2])
				);
    	}
}

double Cfermi_distribution::CalcEquilibriumDensity_sp(double T, double mu) {

    /*  if (T > 0.) {

          auto integrand = [ = ] (double Energy) ->double {

              return GetDOS_sp(Energy) * GetFermiDistribution(Energy, T, mu);
          };


          return integrator_ptr->integrate_qag(integrand, LowerBound(), UpperBound());

      } else {

          std::cerr << "\nPlease start calculation at temperature greater than zero!\n";

          std::exit(1);
      }*/
    	double sum = 0.0;

    	if (T > 0.0) {
	
        	for (unsigned int j = 1; j < n - 1; j++) {

            	sum += DOS_sp[j] * GetFermiDistribution(energy[j], T, mu) * (energy[j + 1] - energy[j - 1]);
        	}
        	
		return 0.5 * (sum
                	+ DOS_sp[0] * GetFermiDistribution(energy[0], T, mu)*(energy[1] - energy[0])
                	+ DOS_sp[n - 1] * GetFermiDistribution(energy[n - 1], T, mu) * (energy[n - 1] - energy[n - 2])
                	);
    	} else {

        	for (unsigned int j = 1; j < n - 1; j++) {
            
			sum += (energy[j] > fermi_energy) ? 0.0 : DOS_sp[j]*(energy[j + 1] - energy[j - 1]);
        	}
       
	 	return 0.5 * (sum + DOS_sp[0]*(energy[1] - energy[0]) + DOS[n - 1]*(energy[n - 1] - energy[n - 2]));
    	} 
}

double Cfermi_distribution::CalcEquilibriumDensity_d(double T, double mu) {
/*
      if (T > 0.) {

          auto integrand = [ = ] (double Energy) ->double {

              return GetDOS_d(Energy) * GetFermiDistribution(Energy, T, mu);
          };


          return integrator_ptr->integrate_qag(integrand, LowerBound(), UpperBound());

      } else {

          std::cerr << "\nPlease start calculation at temperature greater than zero!\n";

          std::exit(1);
      }*/

    	double sum = 0.0;

    	if (T > 0.0) {

        	for (unsigned int j = 1; j < n - 1; j++) {

            		sum += DOS_d[j] * GetFermiDistribution(energy[j], T, mu) * (energy[j + 1] - energy[j - 1]);
        	}
        
		return 0.5 * (sum
                	+ DOS_d[0] * GetFermiDistribution(energy[0], T, mu)*(energy[1] - energy[0])
                	+ DOS_d[n - 1] * GetFermiDistribution(energy[n - 1], T, mu) * (energy[n - 1] - energy[n - 2])
                	);
    	} else {

        	for (unsigned int j = 1; j < n - 1; j++) {

            	sum += (energy[j] > fermi_energy) ? 0.0 : DOS_d[j]*(energy[j + 1] - energy[j - 1]);
        	}

        	return 0.5 * (sum + DOS_d[0]*(energy[1] - energy[0]));
    	}
}

/*
double Cfermi_distribution::df_dTe(double energy, double T, double mu) {

    	return ((energy - mu) / (k_B_au * T)) / (2.0 + 2.0 * std::cosh((energy - mu) / (k_B_au * T)));
}

double Cfermi_distribution::df_dmu(double energy, double T, double mu) {

    	return (1.0 / (k_B_au * T)) / (2.0 + 2.0 * std::cosh((energy - mu) / (k_B_au * T)));
}
*/

double Cfermi_distribution::Get_CT(double T, double mu) {

    	double sum = 0.0;

    	if (T == 0.0) {

        	return 0.0;

    	} else {
        	for (unsigned int i = 1; i < n - 1; i++) {
				double sum0 = (energy[i] * DOS[i]*(energy[i] - mu) / (2.0 + 2.0 * cosh((energy[i] - mu) / (k_B_au * T))))*(energy[i + 1] - energy[i - 1]);
				if(std::isfinite(sum0)){
            		sum += sum0;
				}
			}

			double sum0 = (energy[0] * DOS[0]*(energy[0] - mu) / (2.0 + 2.0 * cosh((energy[0] - mu) / (k_B_au * T))))*(energy[1] - energy[0]);
			if(std::isfinite(sum0)){
            	sum += sum0;
			}

			sum0 = (energy[n - 1] * DOS[n - 1]*(energy[n - 1] - mu) / (2.0 + 2.0 * cosh((energy[n - 1] - mu) / (k_B_au * T))))*(energy[n - 1] - energy[n - 2]);
			if(std::isfinite(sum0)){
            	sum += sum0;
			}
			
			return 0.5 / (k_B_au * T * T)*(sum);
    	}
    	
}

double Cfermi_distribution::Get_CT_sp(double T, double mu) {


    	double sum = 0.0;

    	if (T == 0.0) {

        	return 0.0;

    	} else {

        	for (unsigned int i = 1; i < n - 1; i++) {
				double sum0 = (energy[i] * DOS_sp[i]*(energy[i] - mu) / (2.0 + 2.0 * cosh((energy[i] - mu) / (k_B_au * T))))*(energy[i + 1] - energy[i - 1]);
				if(std::isfinite(sum0)){
            		sum += sum0;
				}
			}

			double sum0 = (energy[0] * DOS_sp[0]*(energy[0] - mu) / (2.0 + 2.0 * cosh((energy[0] - mu) / (k_B_au * T))))*(energy[1] - energy[0]);
			if(std::isfinite(sum0)){
            	sum += sum0;
			}

			sum0 = (energy[n - 1] * DOS_sp[n - 1]*(energy[n - 1] - mu) / (2.0 + 2.0 * cosh((energy[n - 1] - mu) / (k_B_au * T))))*(energy[n - 1] - energy[n - 2]);
			if(std::isfinite(sum0)){
            	sum += sum0;
			}

			return 0.5 / (k_B_au * T * T)*(sum);
    	}
}

double Cfermi_distribution::Get_CT_d(double T, double mu) {

    	double sum = 0.0;

    	if (T == 0.0) {

        	return 0.0;

    	} else {

        	for (unsigned int i = 1; i < n - 1; i++) {
				double sum0 = ((energy[i] * DOS_d[i]*(energy[i] - mu) / (2.0 + 2.0 * cosh((energy[i] - mu) / (k_B_au * T))))*(energy[i + 1] - energy[i - 1]));
				if(std::isfinite(sum0)){
            		sum += sum0;
				}
        	}

			double sum0 = (energy[0] * DOS_d[0]*(energy[0] - mu) / (2.0 + 2.0 * cosh((energy[0] - mu) / (k_B_au * T))))*(energy[1] - energy[0]);
			if(std::isfinite(sum0)){
            	sum += sum0;
			}

			sum0 = (energy[n - 1] * DOS_d[n - 1]*(energy[n - 1] - mu) / (2.0 + 2.0 * cosh((energy[n - 1] - mu) / (k_B_au * T))))*(energy[n - 1] - energy[n - 2]);
			if(std::isfinite(sum0)){
            	sum += sum0;
			}

        	return 0.5 / (k_B_au * T * T)*(sum);
    	}
}

double Cfermi_distribution::Get_CT_mu(double T, double mu) {

    	double sum = 0.0;

    	if (T == 0.0) {

        	return 0.0;

    	} else {

        	for (unsigned int i = 1; i < n - 1; i++) {
				double sum0 = ((energy[i] * DOS[i] / (2.0 + 2.0 * cosh((energy[i] - mu) / (k_B_au * T))))*(energy[i + 1] - energy[i - 1]));
				if(std::isfinite(sum0)){
            		sum += sum0;
				}
			}

			double sum0 = (energy[0] * DOS[0] / (2.0 + 2.0 * cosh((energy[0] - mu) / (k_B_au * T))))*(energy[1] - energy[0]);
			if(std::isfinite(sum0)){
            	sum += sum0;
			}

			sum0 = (energy[n - 1] * DOS[n - 1] / (2.0 + 2.0 * cosh((energy[n - 1] - mu) / (k_B_au * T))))*(energy[n - 1] - energy[n - 2]);
			if(std::isfinite(sum0)){
            	sum += sum0;
			}

			return 0.5 / (k_B_au * T)*(sum);
    	}
}

double Cfermi_distribution::Get_CT_mu_sp(double T, double mu) {

  
    	double sum = 0.0;

    	if (T == 0.0) {

        	return 0.0;

    	} else {

        	for (unsigned int i = 1; i < n - 1; i++) {
				double sum0 = ((energy[i] * DOS_sp[i] / (2.0 + 2.0 * cosh((energy[i] - mu) / (k_B_au * T))))*(energy[i + 1] - energy[i - 1]));
				if(std::isfinite(sum0)){
            		sum += sum0;
				}
			}

			double sum0 = (energy[0] * DOS_sp[0] / (2.0 + 2.0 * cosh((energy[0] - mu) / (k_B_au * T))))*(energy[1] - energy[0]);
			if(std::isfinite(sum0)){
            	sum += sum0;
			}

			sum0 = (energy[n - 1] * DOS_sp[n - 1] / (2.0 + 2.0 * cosh((energy[n - 1] - mu) / (k_B_au * T))))*(energy[n - 1] - energy[n - 2]);
			if(std::isfinite(sum0)){
            	sum += sum0;
			}

			return 0.5 / (k_B_au * T)*(sum);
    	}
}

double Cfermi_distribution::Get_CT_mu_d(double T, double mu) {

    	double sum = 0.0;

    	if (T == 0.0) {

        	return 0.0;

    	} else {

        	for (unsigned int i = 1; i < n - 1; i++) {
				double sum0 = ((energy[i] * DOS_d[i] / (2.0 + 2.0 * cosh((energy[i] - mu) / (k_B_au * T))))*(energy[i + 1] - energy[i - 1]));
				if(std::isfinite(sum0)){
            		sum += sum0;
				}
			}

			double sum0 = (energy[0] * DOS_d[0] / (2.0 + 2.0 * cosh((energy[0] - mu) / (k_B_au * T)))) *(energy[1] - energy[0]);
			if(std::isfinite(sum0)){
            	sum += sum0;
			}

			sum0 = (energy[n - 1] * DOS_d[n - 1] / (2.0 + 2 * cosh((energy[n - 1] - mu) / (k_B_au * T))))*(energy[n - 1] - energy[n - 2]);
			if(std::isfinite(sum0)){
            	sum += sum0;
			}

			return 0.5 / (k_B_au * T)*(sum);
    	}
}

double Cfermi_distribution::Get_P(double T, double mu) {

    double sum = 0.0;

	if(T==0.){
		
		return 0.0;
	}else{

		for (unsigned int i = 1; i < n - 1; i++) {
			double sum0 = ((DOS[i]*(energy[i] - mu) / (2.0 + 2.0 * cosh((energy[i] - mu) / (k_B_au * T))))*(energy[i + 1] - energy[i - 1]));
			if(std::isfinite(sum0)){
            	sum += sum0;
			}
		}

		double sum0 = (DOS[0]*(energy[0] - mu) / (2.0 + 2.0 * cosh((energy[0] - mu) / (k_B_au * T)))) *(energy[1] - energy[0]);
		if(std::isfinite(sum0)){
            sum += sum0;
		}

		sum0 = (DOS[n - 1]*(energy[n - 1] - mu) / (2.0 + 2.0 * cosh((energy[n - 1] - mu) / (k_B_au * T))))*(energy[n - 1] - energy[n - 2]);
		if(std::isfinite(sum0)){
            sum += sum0;
		}

		return 0.5 / (k_B_au * T * T)*(sum);
	}
}

double Cfermi_distribution::Get_P_sp(double T, double mu) {

    double sum = 0.0;

	if (T==0.){
		
		return 0.0;
	
	}else{

    	for (unsigned int i = 1; i < n - 1; i++) {
			double sum0 = ((DOS_sp[i]*(energy[i] - mu) / (2.0 + 2.0 * cosh((energy[i] - mu) / (k_B_au * T))))*(energy[i + 1] - energy[i - 1]));
			if(std::isfinite(sum0)){
            	sum += sum0;
			}
		}

		double sum0 = (DOS_sp[0]*(energy[0] - mu) / (2.0 + 2.0 * cosh((energy[0] - mu) / (k_B_au * T))))*(energy[1] - energy[0]);
		if(std::isfinite(sum0)){
            sum += sum0;
		}

		sum0 = (DOS_sp[n - 1]*(energy[n - 1] - mu) / (2.0 + 2.0 * cosh((energy[n - 1] - mu) / (k_B_au * T))))*(energy[n - 1] - energy[n - 2]);
		if(std::isfinite(sum0)){
            sum += sum0;
		}

		return 0.5 / (k_B_au * T * T)*(sum);

	}
}

double Cfermi_distribution::Get_P_d(double T, double mu) {

    double sum = 0.;

	if (T==0.){
		
		return 0.;

	}else{

    	for (unsigned int i = 1; i < n - 1; i++) {
			double sum0 = ((DOS_d[i]*(energy[i] - mu) / (2.0 + 2.0 * cosh((energy[i] - mu) / (k_B_au * T))))*(energy[i + 1] - energy[i - 1]));
			if(std::isfinite(sum0)){
            	sum += sum0;
			}
		}

		double sum0 = (DOS_d[0]*(energy[0] - mu) / (2.0 + 2.0 * cosh((energy[0] - mu) / (k_B_au * T)))) *(energy[1] - energy[0]);
		if(std::isfinite(sum0)){
            sum += sum0;
		}

		sum0 = (DOS_d[n - 1]*(energy[n - 1] - mu) / (2.0 + 2.0 * cosh((energy[n - 1] - mu) / (k_B_au * T))))*(energy[n - 1] - energy[n - 2]);
		if(std::isfinite(sum0)){
            sum += sum0;
		}

		return 0.5 / (k_B_au * T * T)*(sum);

	}
}

double Cfermi_distribution::Get_P_mu(double T, double mu) {

   double sum = 0.;

	if (T==0){
	
		return 0.0;

	}else{

    	for (unsigned int i = 1; i < n - 1; i++) {
			double sum0 = ((DOS[i] / (2.0 + 2.0 * cosh((energy[i] - mu) / (k_B_au * T))))*(energy[i + 1] - energy[i - 1]));
			if(std::isfinite(sum0)){
            	sum += sum0;
			}
		}

		double sum0 = (DOS[0] / (2.0 + 2.0 * cosh((energy[0] - mu) / (k_B_au * T))))*(energy[1] - energy[0]);
		if(std::isfinite(sum0)){
            sum += sum0;
		}

		sum0 = (DOS[n - 1] / (2.0 + 2 * cosh((energy[n - 1] - mu) / (k_B_au * T))))*(energy[n - 1] - energy[n - 2]);
		if(std::isfinite(sum0)){
            sum += sum0;
		}

		return 0.5 / (k_B_au * T)*(sum);

	}
}

double Cfermi_distribution::Get_P_mu_sp(double T, double mu) {

    double sum = 0.;

	if (T==0){
	
		return 0.0;

	}else{

    	for (unsigned int i = 1; i < n - 1; i++) {
			double sum0 = ((DOS_sp[i] / (2.0 + 2.0 * cosh((energy[i] - mu) / (k_B_au * T))))*(energy[i + 1] - energy[i - 1]));
			if(std::isfinite(sum0)){
            	sum += sum0;
			}
		}

		double sum0 = (DOS_sp[0] / (2.0 + 2.0 * cosh((energy[0] - mu) / (k_B_au * T))))*(energy[1] - energy[0]);
		if(std::isfinite(sum0)){
            sum += sum0;
		}

		sum0 = (DOS_sp[n - 1] / (2.0 + 2.0 * cosh((energy[n - 1] - mu) / (k_B_au * T))))*(energy[n - 1] - energy[n - 2]);
		if(std::isfinite(sum0)){
            sum += sum0;
		}

		return 0.5 / (k_B_au * T)*(sum);

	}
}

double Cfermi_distribution::Get_P_mu_d(double T, double mu) {

    double sum = 0.;

	if (T==0){
	
		return 0.;

	}else{

    	for (unsigned int i = 1; i < n - 1; i++) {
			double sum0 = ((DOS_d[i] / (2.0 + 2.0 * cosh((energy[i] - mu) / (k_B_au * T))))*(energy[i + 1] - energy[i - 1]));
			if(std::isfinite(sum0)){
            	sum += sum0;
			}
		}

		double sum0 = (DOS_d[0] / (2.0 + 2.0 * cosh((energy[0] - mu) / (k_B_au * T))))*(energy[1] - energy[0]);
		if(std::isfinite(sum0)){
			sum += sum0;
		}

		sum0 = (DOS_d[n - 1] / (2.0 + 2 * cosh((energy[n - 1] - mu) / (k_B_au * T))))*(energy[n - 1] - energy[n - 2]);
		if(std::isfinite(sum0)){
            sum += sum0;
		}

		return 0.5 / (k_B_au * T)*(sum);
	
	}
}

double Cfermi_distribution::heatCapacity(double T, double mu){

	return Get_CT(T, mu) - ((Get_CT_mu(T, mu)*Get_P(T, mu))/Get_P_mu(T, mu));
}


double Cfermi_distribution::heatCapacity_sp(double T, double mu){

	return Get_CT_sp(T, mu) - ((Get_CT_mu_sp(T, mu)*Get_P_sp(T, mu))/Get_P_mu_sp(T, mu));
}


double Cfermi_distribution::heatCapacity_d(double T, double mu){

	return Get_CT_d(T, mu) - ((Get_CT_mu_d(T, mu)*Get_P_d(T, mu))/Get_P_mu_d(T, mu));
}



















