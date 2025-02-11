#ifndef ROOT_FINDING_H_
#define ROOT_FINDING_H_

#include<iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

#include"gsl_function_pp.h"
#include"gsl_multiroot_function_pp.h"
#include"physical_constants.h"

#include <gsl/gsl_vector.h>
#include <stdlib.h>
#include <stdio.h>

double onedimensional_root_finding(const auto& func, double x_lo, double x_hi) {
    gsl_function_pp<decltype(func) > Fp(func);
    gsl_function *F = static_cast<gsl_function*> (&Fp);

    int status;
    int iter = 0, max_iter = 10000;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r = 0;


    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, F, x_lo, x_hi);

    do {
        iter++;
        //if(iter < 1)
        //std::cout <<  iter << " " <<  s->x_lower << " " << s->x_upper <<std::endl;
                    

        status = gsl_root_fsolver_iterate(s);
        r = gsl_root_fsolver_root(s);
        x_lo = gsl_root_fsolver_x_lower(s);
        x_hi = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(x_lo, x_hi,
                0, 1e-13);

    } while (status == GSL_CONTINUE && iter < max_iter);
    if (iter == max_iter) {
        //EXIT("Maximum iterations reached");
    }
    gsl_root_fsolver_free(s);

    return r;
}

double multiroot_finding(const auto& func, double* xinit) {
    
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
    
    int status;
    size_t iter = 0;
    const size_t dimensions = 2;
    gsl_vector *x = gsl_vector_alloc(dimensions);
    gsl_multiroot_function_pp<decltype(func) > Fp(func, dimensions);
    gsl_multiroot_function *F = static_cast<gsl_multiroot_function*> (&Fp);

    gsl_vector_set(x, 0, xinit[0]);
    gsl_vector_set(x, 1, xinit[1]);
    
    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T, 2);
    

        /*printf ("iter = %3u x = % .3f % .3f "
              "f(x) = % .3e % .3e\n",
              iter,
              gsl_vector_get (s->x, 0), 
              gsl_vector_get (s->x, 1),
              gsl_vector_get (s->f, 0), 
              gsl_vector_get (s->f, 1));*/
    gsl_multiroot_fsolver_set(s, F, x);
    do {
        iter++;
        //        printf ("iter = %3u x = % .3f % .3f "
        //          "f(x) = % .3e % .3e\n",
        //          iter,
        //          gsl_vector_get (s->x, 0), 
        //          gsl_vector_get (s->x, 1),
        //          gsl_vector_get (s->f, 0), 
        //          gsl_vector_get (s->f, 1));

        status = gsl_multiroot_fsolver_iterate(s);

        //        printf ("iter = %3u x = % .3f % .3f "
        //          "f(x) = % .3e % .3e\n",
        //          iter,
        //          gsl_vector_get (s->x, 0), 
        //          gsl_vector_get (s->x, 1),
        //          gsl_vector_get (s->f, 0), 
        //          gsl_vector_get (s->f, 1));

        if (status) //
            break;

        status =
                gsl_multiroot_test_residual(s->f, 1e-9);
    } while (status == GSL_CONTINUE && iter < 1000);
    if (iter >= 1000) {
        std::cerr << "ERROR: Maximum number of iteration exceeded." << std::endl;
    }
    xinit[0] = gsl_vector_get(s->x, 0);
    xinit[1] = gsl_vector_get(s->x, 1);
    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);
    return 0;

}

#endif
