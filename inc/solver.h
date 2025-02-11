/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   solver.h
 * Author: Ndione, Ag Rethfeld
 *
 * Created on 12. Juni 2018, 17:33
 */

#ifndef SOLVER_H
#define SOLVER_H

#include<iostream>
#include <fstream>
#include<boost/array.hpp>
#include<boost/numeric/odeint.hpp>
#include<boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/filesystem.hpp>

#include "ode.h"
#include "dielectric_function.h"
#include "Laser.h"
#include "fermi_distribution.h"
#include "atom_unit_physics.h"
#include "physical_constants.h"
#include "TwoBandWithPhonons.h"
#include "TwoBandWithPhonons2T.hpp"
#include "ThreeBandWithPhonons.h"
#include "ThreeBandWithPhonons2T.hpp"
#include "FourBandWithPhonons.h"
#include "FourBandWithPhonons2T.hpp"
#include "used_data.h"



using namespace boost::filesystem;
using namespace boost::numeric::odeint;


typedef std::vector<double> dv;

class solver {

	public:
		solver();
    		solver(int index);
    		virtual  ~solver();

    		void run();
   		void odeintsolver_seb(const dv &q, dv &q_dot, double time);

    		void Print(const std::vector<double> &x, const double time);
   		void WritetoFile(const std::vector<double> &x, const double time, int);
		int getChosenModel(); 



	private:

		ODE* ode_point;
    		std::ofstream outputfile;
		std::ofstream inputfile;
		std::ofstream paramFile;
   		Cfermi_distribution *fermi_ptr;
   		Claser *laser_point;
		int chooseModel;

};



#endif /* SOLVER_H */

