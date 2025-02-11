/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   qag_integration.h
 * Author: ndione
 *
 * Created on March 2, 2017, 11:37 AM
 */

#ifndef QAG_INTEGRATION_H
#define QAG_INTEGRATION_H


class Cintegrate{
    
public:
	Cintegrate()
	{
		integrate_workspace = gsl_integration_workspace_alloc (11000);
	}
	~Cintegrate(){
		gsl_integration_workspace_free(integrate_workspace);
	}
	
	double integrate(const auto& func, double min, double max)
	{
		gsl_function_pp<decltype(func)> Fp(func);
		gsl_function *F = static_cast<gsl_function*>(&Fp);

        	double estimated_error=-1;
        	double integration_result=NAN;
        	int intervals=11000;

        	double eps_abs=0;
        	double eps_rel=1e-8;
        	int integral_type=GSL_INTEG_GAUSS61;

        	error=gsl_integration_qag(F,min,max,eps_abs,eps_rel,intervals,6,integrate_workspace,&integration_result,&estimated_error);
        	return integration_result;
	}

private:
	gsl_integration_workspace* integrate_workspace;
	int error;
	

};

#endif /* QAG_INTEGRATION_H */

