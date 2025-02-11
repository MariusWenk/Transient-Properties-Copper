#ifndef INTEGRATION_H_
#define INTEGRATION_H_
#include"gsl_function_pp.h"
#include<gsl/gsl_integration.h>
#include "physical_constants.h"


class Cintegration {
public:

    //!Constructor
    Cintegration();
    //!Destructor
    ~Cintegration();

    /*!********************************************************************
     * \brief Base integration function using adaptive Gauss-Kronrod quadrature routine
     *
     * Hands over the function to integrate to integrate_qag.
     * The function has to be defined here, due to the auto in the signature.
     * 
     * \param func	Address of the lambda function, that should be integrated
     * \param min	Lower border of the integration interval
     * \param max	higher border of the integration interval
     * \return Integrated value of Cintegration::integrate_qag
     * 
     **********************************************************************/
    double integrate(const auto& func, double min, double max) {
        //return simpson(func,min,max,10000);
        return integrate_qag(func, min, max);
    }

    /*!********************************************************************
     * \brief Integration function using adaptive Gauss-Kronrod quadrature routine
     *
     * The function has to be defined here, due to the auto in the signature.
     *
     * \param func	Address of the lambda function, that should be integrated
     * \param min	Lower border of the integration interval
     * \param max	higher border of the integration interval
     * \return Integrated
     *
     **********************************************************************/
    double integrate_qag(const auto& func, double min, double max) {
        gsl_function_pp<decltype(func) > Fp(func);
        gsl_function *F = static_cast<gsl_function*> (&Fp);

        double integration_result = NAN;

        last_error = gsl_integration_qag(F, min, max, eps_abs, eps_rel, interval_max, qag_integral_type, integrate_workspace, &integration_result, &estimated_error);
        return integration_result;

    }

    /*!********************************************************************
     * \brief Integration function using Cauchys principle value
     *
     * This function can integrate functions with one known, (easy) singularity.
     * Therefore the integrand \f$f(x)\f$ must be splitted and the singularity
     * \f$1/(x-c)\f$ must be seperated, that \f$f(x) = func * 1/(x-c)\f$
     *
     * This function has to be defined here, due to the auto in the signature.
     *
     * \param func	Address of the lambda function, that should be integrated
     * \param min	Lower border of the integration interval
     * \param max	higher border of the integration interval
     * \return Cauchys principle value
     **********************************************************************/
    double integrate_qawc(const auto& func, double min, double max, double singularity) {

        gsl_function_pp<decltype(func) > Fp(func);
        gsl_function *F = static_cast<gsl_function*> (&Fp);

        double integration_result = NAN;

        last_error = gsl_integration_qawc(F, min, max, singularity, eps_abs, eps_rel, interval_max, integrate_workspace, &integration_result, &estimated_error);
        return integration_result;

    }
    double simpson(const auto& y_func, double start, double ende, int number_of_disc_points) {
    double result = 0.0;
    double x = (ende - start) / number_of_disc_points;
    for (int i = 0; i != number_of_disc_points; i++) {
        double c = ((x * i) + (x * (i + 1))) / 2.0;
        double h3 = fabs(x * (i + 1) - (x * i)) / 6.0;
        result += h3 * (y_func(x * i) + 4.0 * y_func(c) + y_func(x * (i + 1)));

    }
    return result;
}
private:
    //!Workspace of the gsl_integration routines
    gsl_integration_workspace* integrate_workspace;

    //!Size of the gsl_integration_workspace
    int interval_max;
    //!Error-value of the GSL-integration
    int last_error;
    //!Estimated error of the GSL-integration
    double estimated_error;
    //!Routine that is chosen in the QAG-Integration
    int qag_integral_type;

    //!Maximum absolute error
    double eps_abs;
    //!Maximum relative error
    double eps_rel;


};
#endif