#include"integration.h"

Cintegration::Cintegration()
{

       	interval_max = 1000;
        estimated_error=-1;
        qag_integral_type=GSL_INTEG_GAUSS61;

        eps_abs=1;  //maximum absolute error
		//eps_abs=0;
        eps_rel=1e-6; //maximum relative error

	integrate_workspace = gsl_integration_workspace_alloc (interval_max);

	
}
Cintegration::~Cintegration()
{
	gsl_integration_workspace_free(integrate_workspace);
}	
//double Cintegration::integrate(const auto& func, double min, double max)
//{
//	//QAG is the default integration
//	return integrate_qag(func, min, max);
//}
//double Cintegration::integrate_qag(const auto& func, double min, double max)
//{
//	gsl_function_pp<decltype(func)> Fp(func);
//	gsl_function *F = static_cast<gsl_function*>(&Fp);
//
//	double integration_result=NAN;
//
//	last_error=gsl_integration_qag(F,min,max,eps_abs,eps_rel,interval_max,qag_integral_type,integrate_workspace,&integration_result,&estimated_error);
//	return integration_result;
//}
