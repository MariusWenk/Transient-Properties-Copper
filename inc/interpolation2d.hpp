#ifndef INTERPOLATION2D_H
#define INTERPOLATION2D_H

//#include<iostream>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "basic_inputparser.hpp"
#include "debug.hpp"

typedef unsigned int uint;

namespace ufd{

    /*! \~english \brief wrapper class for gsl gsl 2D-Interpolation

        This class allows us to do a 2D-Interpolation.
        Pointers, vectors or files can all be used as input.
        Adding the possibility of calculating derivatives is planned for the future.

        \~german \brief Wrapper-Klasse für gsl 2D-Interpolation

        Diese Klasse bietet die Möglichkeit eine 2D-Interpolation durchzuführen.
        Dabei können als Eingaben pointer, Vektoren oder ein File verwendet werden.
        In Zukunft soll die Möglichkeit bestehen sich die Ableitung zurückgeben zu lassen
    */
    class Interpolation2D {
        public:

            //! \~ constructor (does nothing)
            Interpolation2D();

            /*! \~ \brief constructor (pointer-input)

                Create a new Interpolation with pointer input
                @param x_array
                @param y_array
                @param z_array
                @param xarray_length
                @param yarray_length
                The values in the x_array and y_array have to be in ascending order, the values in
                the z_array are accessed like z_ij = z_array[j*xarray_length +i]
                \sa Interp2D(std::vector<double> x_array, std::vector<double> y_array, std::vector<double> z_array)
                \sa Interp2D(std::string filename, uint x_col, uint y_col, uint z_col)
            */
            Interpolation2D(double* x_array, double* y_array, double* z_array, uint xarray_length, uint yarray_length);

            /*! \~ \brief constructor (pointer-input with array of array)

                Create a new Interpolation with pointer input
                @param x_array
                @param y_array
                @param z_array
                @param xarray_length
                @param yarray_length
                The values in the x_array and y_array have to be in ascending order, the values in
                the z_array are accessed like z_ij = z_array[i][j]
                \sa Interp2D(std::vector<double> x_array, std::vector<double> y_array, std::vector<double> z_array)
                \sa Interp2D(std::string filename, uint x_col, uint y_col, uint z_col)
            */
            Interpolation2D(double* x_array, double* y_array, double** z_array, uint xarray_length, uint yarray_length);

            /*! \~ \brief constructor (vector-input)

                Create a new Interpolation with vector input
                @param x_array
                @param y_array
                @param z_array
                The values in the x_array and y_array have to be in ascending order, the values in
                the z_array are accessed like z_ij = z_array[j*xarray_length +i]
                \sa Interp2D(double* x_array, double* y_array, double* z_array, uint xarray_length, uint yarray_length)
                \sa Interp2D(std::string filename, uint x_col, uint y_col, uint z_col)
            */
            Interpolation2D(std::vector<double> x_array, std::vector<double> y_array, std::vector<double> z_array);

            /*! \~ \brief constructor (vector-input with vector of vectors)

                Create a new Interpolation with vector input
                @param x_array
                @param y_array
                @param z_array
                The values in the x_array and y_array have to be in ascending order, the values in
                the z_array are accessed like z_ij = z_array[i][j]
                \sa Interp2D(double* x_array, double* y_array, double* z_array, uint xarray_length, uint yarray_length)
                \sa Interp2D(std::string filename, uint x_col, uint y_col, uint z_col)
                \sa Interp2D(std::vector<double> x_array, std::vector<double> y_array, std::vector<double> z_array)
            */
            Interpolation2D(std::vector<double> x_array, std::vector<double> y_array, std::vector<std::vector<double>> z_array);

            /*! \~ \brief Konstruktor (file-input)

                Create a new Interpolation with file input
                the function reads the data from \a filename with the data for x,y,z in the specified columns
                and creates an interpolation
                the file should have the following format (seperated by tabs):
                x_1  y_1  z_11
                x_1  y_2  z_12
                ...
                x_1  y_n  z_1n
                x_2  y_1  z_21
                x_2  y_2  z_22
                ...
                x_m  y_n  z_mn
                @param filename
                @param x_col
                @param y_col
                @param z_col
                \sa Interp2D(double* x_array, double* y_array, double* z_array, uint xarray_length, uint yarray_length)
                \sa Interp2D(std::vector<double> x_array, std::vector<double> y_array, std::vector<double> z_array)
            */
            Interpolation2D(std::string filename, uint x_col, uint y_col, uint z_col);

            //! \~ destructor
            ~Interpolation2D();

            /*! \~english \brief returns interpolated value
                \~german \brief Rückgabe eines interpolierten Wertes

                \~
                @param x_value
                @param y_value
                @return f(x,y)
            */
            double getValue(double x_value, double y_value) const;

            /*
             *
             * @param xvalue
             * @return f'(x)
             */
            //    double get_deriv(double xvalue, double yvalue) const; //derivation only in special directions possible, because whole derivation is a matrix
            /*
             *
             * @param xvalue
             * @return f''(x)
             */
            //    double get_deriv2(double xvalue, double yvalue) const; //derivation only in special directions possible, because whole derivation is a matrix

            /*! \~english \brief Update the z Values of this interpolation with pointer input
                \~german \brief Update the z Values of this interpolation with pointer input

                \~
                @param z_array
                \sa update(std::vector<double> z_array)
                \sa update(double** z_array)
                \sa update(std::vector<std::vector<double>> z_array)
            */
            void update(double* z_array);

            /*! \~english \brief Update the z Values of this interpolation with 2d pointer
                \~german \brief Update the z Values of this interpolation with 2d pointer

                \~
                @param z_array
                \sa update(std::vector<double> z_array)
                \sa update(double* z_array)
                \sa update(std::vector<std::vector<double>> z_array)
            */
            void update(double** z_array);

            /*! \~english \brief Update the z Values of this interpolation with vector input
                \~german \brief Update the z Values of this interpolation with vector input

                \~
                @param z_array
                \sa update(double* z_array)
                \sa update(double** z_array)
                \sa update(std::vector<std::vector<double>> z_array)
            */
            void update(std::vector<double> z_array);

            /*! \~english \brief Update the z Values of this interpolation with 2d vector input
                \~german \brief Update the z Values of this interpolation with 2d vector input

                \~
                @param z_array
                \sa update(double* z_array)
                \sa update(double** z_array)
                \sa update(std::vector<double> z_array)
            */
            void update(std::vector<std::vector<double>> z_array);

            /*! \~english \brief returns smallest x-value for the interpolation
                \~german \brief Rückgabe des minimalen x-Wertes für die Interpolation

                \~
                @return minimal x value
            */
            double getMinX() const;

            /*! \~english \brief returns largest x-value for the interpolation
                \~german \brief Rückgabe des maximalen x-Wertes für die Interpolation

                \~
                @return maximal x value
            */
            double getMaxX() const;

            /*! \~english \brief returns smallest y-value for the interpolation
                \~german \brief Rückgabe des minimalen y-Wertes für die Interpolation

                \~
                @return minimal y value
            */
            double getMinY() const;

            /*! \~english \brief returns largest y-value for the interpolation
                \~german \brief Rückgabe des maximalen y-Wertes für die Interpolation

                \~
                @return maximal y value
            */
            double getMaxY() const;

            /*! \~english \brief returns x-input length
                \~german \brief Rückgabe der Länge der X-Inputdaten

                \~
                @return x_array_length_
            */
            uint getXLength() const;

            /*! \~english \brief returns y-input length
                \~german \brief Rückgabe der Länge der Y-Inputdaten

                \~
                @return y_array_length_
            */
            uint getYLength() const;

            /*! \~english \brief returns x-input array
                \~german \brief Rückgabe des Input-x-Arrays

                \~
                @return x_array_
            */
            double* getXArray() const;

            /*! \~english \brief returns y-input array
                \~german \brief Rückgabe des Input-y-Arrays

                \~
                @return y_array_
            */
            double* getYArray() const;

            /*! \~english \brief returns z-input array
                \~german \brief Rückgabe des Input-z-Arrays

                \~
                @return z_array_
            */
            double* getZArray() const;


        private:

            /*! \~english \brief starts interpolation process

                initializes interpolation and calls \a update();
                called by constructors

                \~german \brief startet den Interpolations-Vorgang

                initialisiert die Interpolation und ruft \a update() auf;
                wird von den Konstruktoren aufgerufen
             */
            void createInterpolation(double* z_array);

            //! \~english saves x-values of the interpolation \~german speichert die x-Werte der Interpolation
            double* x_array_;

            //! \~english saves y-values of the interpolation \~german speichert die y-Werte der Interpolation
            double* y_array_;

            //! \~english saves z-values of the interpolation \~german speichert die z-Werte der Interpolation
            double* z_array_;

            //! \~english saves length of x-input \~german speichert die Länge der X-Inputdaten
            uint x_array_length_;
            //! \~english saves length of y-input \~german speichert die Länge der Y-Inputdaten
            uint y_array_length_;

            // \~ gsl object
            const gsl_interp2d_type *T = gsl_interp2d_bilinear;
            // \~ gsl object
            gsl_interp_accel* xacc_;
            // \~ gsl object
            gsl_interp_accel* yacc_;
            // \~ gsl object
            gsl_spline2d* spline_;
    };
}
#endif /* INTERPOLATION2D_H */
