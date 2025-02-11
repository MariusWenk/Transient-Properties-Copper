/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   integration_routines.cpp
 * Author: ndione
 * 
 * Created on 1. Februar 2018, 12:40
 */

#include "integration_routines.h"

Cintegration_routines::Cintegration_routines() {

}

double Cintegration_routines::trapez(double *x, double *y, int number_of_disc_points) {
    double ergebnis = 0.0;
    for (int i = 0; i != number_of_disc_points - 1; i++) {
        ergebnis += (x[i + 1] - x[i])*(y[i + 1] + y[i]) / 2;
    }
    return ergebnis;

}

double Cintegration_routines::trapez(const auto& x_func, const auto& y_func, int number_of_disc_points) {
    double ergebnis = 0.0;
    for (int i = 0; i != number_of_disc_points; i++) {

        ergebnis += (x_func(i + 1) - x_func(i))*(y_func(i + 1) + y_func(i)) / 2;
    }

    return ergebnis;

}

double Cintegration_routines::simpson(const auto& x_func, const auto& y_func, int number_of_disc_points) {
    double result = 0;
    for (int i = 0; i != number_of_disc_points; i++) {
        double c = (x_func(i) + x_func(i + 1)) / 2.0;
        double h3 = abs(x_func(i + 1) - x_func(i)) / 6.0;
        result += h3 * (y_func(i) + 4.0 * y_func(c) + y_func(i + 1));
    }
    return result;
}

double Cintegration_routines::trapez(double* y_func, double start, double ende, int number_of_disc_points) {
    double ergebnis = 0.0;
    double x = (ende - start) / number_of_disc_points;
    for (int i = 0; i != number_of_disc_points; i++) {

        ergebnis += (x * (i + 1) - (x * i))*(y_func[(i + 1)] + y_func[i]) / 2;

    }

    return ergebnis;

}

double Cintegration_routines::simpson(const auto& y_func, double start, double ende, int number_of_disc_points) {
    double result = 0.0;
    double x = (ende - start) / number_of_disc_points;
    for (int i = 0; i != number_of_disc_points; i++) {
        double c = ((x * i) + (x * (i + 1))) / 2.0;
        double h3 = fabs(x * (i + 1) - (x * i)) / 6.0;
        result += h3 * (y_func(x * i) + 4.0 * y_func(c) + y_func(x * (i + 1)));

    }
    return result;
}

