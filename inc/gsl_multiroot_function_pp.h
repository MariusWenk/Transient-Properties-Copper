#ifndef GSL_MULTIROOT_FUNCTION_PP_H
#define GSL_MULTIROOT_FUNCTION_PP_H

#include <gsl/gsl_multiroots.h> // defines gsl_multiroot_function

/*!**************************************************************************************
	\brief Wrapper class to transform a class function into a gsl_multiroot_function.

	This wrapper class is necessary because gsl_multiroot_functions have to be static
	which causes a problem when desiring to use class functions. <br/>
	Hence, it transforms the provided function into a gsl_multiroot_function consisting 
	of int(* f)(const gsl_vector * x, void * params, gsl_vector * f), size_t n and 
	void * params, where the function f stores its results for argument x and parameters
	params in the gsl_vector f, n is the dimension of the system and params is the
	pointer to the parameters of the function.<br/>
	The constructor saves the provided pointer to a to a constant class variable 
	\ref _func. It then assigns the reference to the private function 
	\ref invoke to gsl_multiroot_function.f, the pointer to this 
	object to the void pointer gsl_multiroot_function.params and the given dimension to 
	gsl_multiroot_function.n. The function \ref invoke in turn returns a static cast of 
	the multiroot_function \ref _func thus making it usable as gsl_multiroot_function.
	<br/>
	To use this wrapper, create a class function, evaluating the functions to search
	root of. This class function needs the type int func_eval(const gsl_vector * x, 
	gsl_vector * f), has to evaluate the functions for the parameters stored in 
	const gsl_vector * x and has to store the result in gsl_vector * f. It then returns
	GSL_SUCCESS. To cast it to a gsl_multiroot_function using this wrapper, we need a
	lambda function. Here, you find an example:<br/>
	auto ptr = [&](const gsl_vector * x, gsl_vector * f)->int{return func_eval(x,f);};
	<br/>
	gsl_multiroot_function_pp<decltype(ptr)> Fp(ptr, dimensions); <br/>
	gsl_multiroot_function * F = static_cast<gsl_multiroot_function * >(&Fp);<br/>
	Thanks to Sebastian Weber for providing \ref gsl_function_pp as basis for this 
	wrapper.
*****************************************************************************************/

template< typename F >
	/* gsl_multiroot_function is declared in gsl/gsl_multiroots.h.
	 * It consists of int(* f)(const gsl_vector x, void *params, gsl_vector f), size_t n and void *params.
	 */
	class gsl_multiroot_function_pp : public gsl_multiroot_function {

	public:
		/*!******************************************************************************
		 * \brief Constructor
		 *
		 * Asigns provided reference to a (class-)function to local constant 
		 * \ref _func. Then assigns reference to private function 
		 * invoke(gsl_vector *x, void* params, gsl_vector f) to gsl_multiroot_function.f, 
		 * the pointer to this object to gsl_multiroot_function.params and the given 
		 * dimension to gsl_multiroot_function.n.
		 * \param func 			reference to (class-)multiroot_function in question;
		 * 						needs the signature int(const gsl_vector, gsl_vector).
		 * \param dimensions	dimensions of the system 
		 ********************************************************************************/
		gsl_multiroot_function_pp(const F& func, size_t dimensions) : _func(func) {
			// Setting _func has to be done using this syntax because _func is const and cannot be set once the class is created.
			// saves pointer to invoke in multiroot_function declared in gsl_multiroot_function
			f = &gsl_multiroot_function_pp::invoke;
			// saves pointer to this class in params
			params = this;
			// saves dimension
			n = dimensions;
		}
	
	private:
		//! Reference to (class-)multiroot_function in question
		const F& _func;
		/*!********************************************************************
		 * \brief Casts function to static making it usable by GSL.
		 *
		 * Casts params->_func(x) which corresponds to this->_func(x) and thus 
		 * the multiroot_function handed over to this class to static which makes it 
		 * usable by the GSL. The parameter list has to match the one of
		 * gsl_multiroot_function.
		 * \param x			argument of multiroot_function in question
		 * \param params	set of parameters
		 * \param f			vector results are saved to
		 * \return static cast of multiroot_function in question
		 **********************************************************************/
		static int invoke(const gsl_vector *x, void *params, gsl_vector *f) {
			// Using params->_func(x)  explicitely is necessary because invoke has to be static for usage with GSL. 
			// Consequently we need to know _func of which instance of the class we have to use.
			return(static_cast<gsl_multiroot_function_pp*>(params)->_func(x, f));
		}
};

#endif // GSL_MULTIROOT_FUNCTION_PP_H
