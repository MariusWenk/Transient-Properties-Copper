/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   integration_routines.h
 * Author: ndione
 *
 * Created on 1. Februar 2018, 12:40
 */

#ifndef INTEGRATION_ROUTINES_H
#define INTEGRATION_ROUTINES_H

#include<cmath>// sqrt,pow
#include<cstdlib>
#include<cstdio>
#include<iostream>

class Cintegration_routines {
public:
    Cintegration_routines();
    /**
     * Trapez-Integration für ein x-Werte-Array und ein y-Werte-Array mit i Stellen 
     * berechnet 0_Integral^i y(x) dx; 
     * @param x_array
     * @param y_array
     * @param number_of_disc_points
     * @return integral
     */
    double trapez(double* x_array, double* y_array, int number_of_disc_points);
    /**
     * Trapez-Integration für eine x-Werte-Funktion und eine y-Werte-Funktion von 0 bis number_of_disc_points
     * berechnet 0_Integral^number_of_disc_points y(x) dx  
     * @param x_func
     * @param y_func
     * @param number_of_disc_points
     * @return integral
     */
    double trapez(const auto& x_func, const auto& y_func, int number_of_disc_points);
    /**
     * Simpson-Integration für eine x-Werte-Funktion und eine y-Werte-Funktion von 0 bis number_of_disc_points
     * berechnet 0_Integral^number_of_disc_points y(x) dx;
     * @param x_func
     * @param y_func
     * @param number_of_disc_points
     * @return integral
     */
    double simpson(const auto& x_func, const auto& y_func, int number_of_disc_points);
    /**
     * Trapez-Integration für eine x-Werte-Funktion und eine y-Werte-Funktion von start bis ende mit number_of_disc_points zwischenstellen
     * berechnet start_Integral^ende y(x) dx ;
     * @param y_func
     * @param start
     * @param ende
     * @param number_of_disc_points
     * @return integral
     */
    double trapez(double* y_func, double start, double ende, int number_of_disc_points);
    /**
     * Trapez-Integration für eine x-Werte-Funktion und eine y-Werte-Funktion von start bis ende mit number_of_disc_points zwischenstellen
     * berechnet start_Integral^ende y(x) dx;
     * @param y_func
     * @param start
     * @param ende
     * @param number_of_disc_points
     * @return integral
     */
    double simpson(const auto& y_func, double start, double ende, int number_of_disc_points);

};

#endif /* INTEGRATION_ROUTINES_H */

