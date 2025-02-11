/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ode.cpp
 * Author: weber
 * 
 * Created on 12. Juni 2018, 17:33
 */

#include "ode.h"

ODE::ODE(){
	
}
ODE::ODE(int number_of_equations)
	:x(number_of_equations), dxdt(number_of_equations){
	is_calculated = false;
}
ODE::~ODE(){

}
dv& ODE::get_change_vector(double time){
	if(!is_calculated){
		calculate_dxdt(time);
		is_calculated = true;
	}
	return dxdt;
}
void ODE::set_new_variables(const dv& changed_variables){
	for (unsigned int i=0;i<changed_variables.size();i++){
		x[i] = changed_variables[i];
	}
	is_calculated=false;
}
void ODE::calculate_dxdt(double time){
	std::cerr << "bla\n";
}

double ODE::get_t_start(){
	return -1;
}

double ODE::get_t_end(){
	return -1;
}

double ODE::get_t_step(){
	return -1;
}

std::string ODE::get_savefilename(){

    return savefilename;
}
