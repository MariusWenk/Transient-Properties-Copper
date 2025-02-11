/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   integrate_template.h
 * Author: ndione
 *
 * Created on March 1, 2017, 1:22 PM
 */

#ifndef INTEGRATE_TEMPLATE_H
#define INTEGRATE_TEMPLATE_H



#ifndef INTEGRATE_H
#define INTEGRATE_H
#include <cmath>
#include <vector>
#include <gsl/gsl_spline.h>
#include <iostream>

template <typename T>
T trapez(T* x, T* y, int length)
{
	T result=0.0;
	for (int i=0;i<length-1;i++){
		result+=(x[i+1]-x[i])*(y[i]+y[i+1]);
	}
	return result*0.5;
}
template <typename T>
T trapez(T* x, T* y, int min, int max){
	T result=0.0;
	for (int i=min;i<max-1;i++){
		result+=(x[i+1]-x[i])*(y[i]+y[i+1]);
	}
	return result*0.5;
}
template <typename T>
T simpson(T* x, T* y, int length){
	gsl_interp_accel* acc;
	gsl_spline* spline;

	acc = gsl_interp_accel_alloc();
	spline = gsl_spline_alloc(gsl_interp_akima,length);
	gsl_spline_init (spline,x,y,length);

	T result=0.0;
	for (int i=0;i<length-1;i++){
		result+=(x[i+1]-x[i])*(y[i]+4*gsl_spline_eval(spline,(x[i+1]+x[i])/2,acc)+y[i+1]);
	}
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
	return result/6;
}

template <typename T>
T pulcherrima(T* x, T* y, int length){
	gsl_interp_accel* acc;
	gsl_spline* spline;

	acc = gsl_interp_accel_alloc();
	spline = gsl_spline_alloc(gsl_interp_akima,length);
	gsl_spline_init (spline,x,y,length);

	T result=0.0;
	for (int i=0;i<length-1;i++){
		T x1=(2*x[i]+x[i+1])/3;
		T x2=(x[i]+2*x[i+1])/3;
		T y1=gsl_spline_eval(spline,x1,acc);
		T y2=gsl_spline_eval(spline,x2,acc);
		result+=(x[i+1]-x[i])*(y[i]+3*y1+3*y2+y[i+1]);
	}
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
	return result/8;
}

template <typename T>
T boole(T* x, T* y, int length){
	gsl_interp_accel* acc;
	gsl_spline* spline;

	acc = gsl_interp_accel_alloc();
	spline = gsl_spline_alloc(gsl_interp_akima,length);
	gsl_spline_init (spline,x,y,length);

	T result=0.0;
	for (int i=0;i<length-1;i++){
		T x1=(3*x[i]+x[i+1])/4;
		T x2=(2*x[i]+2*x[i+1])/4;
		T x3=(x[i]+3*x[i+1])/4;
		T y1=gsl_spline_eval(spline,x1,acc);
		T y2=gsl_spline_eval(spline,x2,acc);
		T y3=gsl_spline_eval(spline,x3,acc);
		result+=(x[i+1]-x[i])*(7*y[i]+32*y1+12*y2+32*y3+7*y[i+1]);
	}
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
	return result/90;
}

template <class T>
T trapez(std::vector<T>* x, std::vector<T>* y, int length)
{
	T result=0.0;
	for (int i=0;i<length-1;i++){
		result+=(x->at(i+1)- x->at(i))*(y->at(i)+y->at(i+1));
	}
	return result*0.5;
}

#endif /* INTEGRATE_H */



#endif /* INTEGRATE_TEMPLATE_H */

