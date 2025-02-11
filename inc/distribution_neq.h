/************************************ 
 * File:   fermi_distribution.h
 * Author: Ndione, Ag Rethfeld
 ************************************
 * Created on 09. January 2019, 10:55
 ************************************/

#ifndef DISTRIBUTION_NEQ_H
#define DISTRIBUTION_NEQ_H

 
#include "fermi_distribution.h"
#include "BoostRootFinding1D.h"
#include "integration.h"

class Cdistribution_neq {

	public:
		Cdistribution_neq();
		Cdistribution_neq(BoostRootFinding1D *boost_pointer, std::string dos_file, double T_elec, double mu_sp, double mu_d);
	        virtual ~Cdistribution_neq();
 
		//!Hands over a fermi distribution for sp-electrons
		double GetDistribution_sp(double);
		//!Hands over a Phi-function for sp-electrons
		double GetPhiFunction_sp(double);
		//!Hands over a fermi distribution for d-electrons
		double GetDistribution_d(double);
		//!Hands over a Phi-function for d-electrons
		double GetPhiFunction_d(double);
		//!Hands over a nonequilibrium distribution derived from two fermi distributions
		//Uses interpolation
		double GetDistributionNonEq(double energy);
		//!Non eq. distribution function using discrete points of DOS
		void GetDistributionNonEq(std::string);
		//!Fermi eq. distribution using discrete points of DOS
		void eq_distribution(std::string);
		//!band-resolved Fermi distribution using discrete points of DOS
		void partialDistribution(std::string);
		//!Hands over nonequilibrium Phi-function
		double GetPhiFunctionNonEq(double energy); 
		//!Hands over the product of a fermi distribution and and associated DOS
		double NumberDensity_sp(double energy);
		double NumberDensity_d(double energy);
		//Check the total number of particles calculated with a nonequilibrium distribution function
		double ParticlesFromDistributionNonEq();
		//Hands over the total number of particles calculated with a fermi distribution
		double ParticlesFromFermiDistribution();
		std::vector<std::vector<double>> load;
 
	 
		double GetTemperature();
		double GetMu_sp();
		double GetMu_d();
	private:

		double T_elec;
		double mu_sp;	
		double mu_d;
 	
		BoostRootFinding1D *boost_pointer;
		 
		Cintegration *integration_pointer;

		double* E;
    		//Array for total DOS
    		double* totDOS;
	    	//Array for sp-DOS
    		double *spDOS;
    		//Array for d-DOS
    		double *dDOS; 

};
#endif /* DISTRIBUTION_NEQ_H */
