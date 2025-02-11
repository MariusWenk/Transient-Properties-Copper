#ifndef GSL_FUNCTION_PP_H
#define GSL_FUNCTION_PP_H

#include <gsl/gsl_math.h> // defines gsl_function

/*!****************************************************************************
	\brief Wrapper class to transform a class function into a gsl_function.

	This wrapper class is necessary because gsl_functions have to be static
	which causes a problem when desiring to use class functions. <br/>
	Hence, it transforms the provided function into a gsl_function consisting of
	double(* function)(double x, void * params) and void * params.<br/>
	The constructor saves the provided pointer to a function to a constant 
	class variable \ref _func. It then assigns the reference to the private 
	function \ref invoke(double x, void * params) to gsl_function.function and 
	the pointer to this object to the void pointer gsl_function.params. The 
	function invoke in turn returns a static cast of the function \ref _func 
	thus making it usable as gsl_function.
	<br/>
	Thanks to Sebastian Weber for providing this wrapper.
*******************************************************************************/

template< typename F >
	/* gsl_function is declared in gsl/gsl_math.h.
	 * It consists of double(* function)(double x, void *params) and void *params.
	 */
	class gsl_function_pp : public gsl_function {

	public:
		/*!********************************************************************
		 * \brief Constructor
		 *
		 * Asigns provided reference to a (class-)function to local constant 
		 * _func. Then assigns reference to private function 
		 * invoke(double x, void* params) to gsl_function.function and the
		 * pointer to this object to gsl_function.params.
		 * \param func reference to (class-)function in question
		 **********************************************************************/
		gsl_function_pp(const F& func) : _func(func) {
			// Setting _func has to be done using this syntax because _func is const and cannot be set once the class is created.
			// saves pointer to invoke in function (hopefully declared in gsl_function)
			function = &gsl_function_pp::invoke;
			// saves pointer to this class in params
			params = this;
		}
	
	private:
		//! Reference to (class-)function in question
		const F& _func;
		/*!********************************************************************
		 * \brief Casts function to static making it usable by GSL.
		 *
		 * Casts params->_func(x) which corresponds to this->_func(x) and thus 
		 * the function handed over to this class to static which makes it 
		 * usable by the GSL. The parameter list has to match the one of
		 * gsl_function.
		 * \param x			argument of function in question
		 * \param params	set of parameters
		 * \return static cast of function in question
		 **********************************************************************/
		static double invoke(double x, void *params) {
			// Using params->_func(x)  explicitely is necessary because invoke has to be static for usage with GSL. 
			// Consequently we need to know _func of which instance of the class we have to use.
			return(static_cast<gsl_function_pp*>(params)->_func(x));
		}
};

#endif // GSL_FUNCTION_PP_H
