/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ode.h
 * Author: weber
 *
 * Created on 12. Juni 2018, 17:33
 */

#ifndef ODE_H
#define ODE_H

/**************************************
 * General class to handle different 
 * first order differential equations
 * of the type dxdt = ....
 **************************************/


#include<iostream>
//#include<vector>
//#include<math.h>
//#include<boost/array.hpp>
#include<boost/numeric/odeint.hpp>
#include "atom_unit_physics.h"
#include "physical_constants.h"
#include "used_data.h"


typedef std::vector<double> dv;

class ODE {

	public:
    		ODE();
    		ODE(int number_of_equations);
    		virtual ~ODE();

    		dv& get_change_vector(double time);
    		void set_new_variables(const dv& changed_variables);
    		virtual double get_t_start();
    		virtual double get_t_end();
    		virtual double get_t_step();
    		std::string get_savefilename();

    		dv x;

	protected:

    		virtual void calculate_dxdt(double time);
		bool is_calculated;
    		std::string savefilename;// = "outputdata.txt";
    		dv dxdt;

};


#endif /* ODE_H */

