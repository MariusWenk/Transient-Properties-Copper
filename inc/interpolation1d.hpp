#ifndef INTERPOLATION1D_H
#define INTERPOLATION1D_H

#include <vector>
#include <gsl/gsl_spline.h>
//#include <iostream>
#include "basic_inputparser.hpp"
#include "debug.hpp"

typedef unsigned int uint;

namespace ufd{

    /*! \~english \brief wrapper class for gsl 1D-Interpolation
    
        This class allows us to do a 1D-Interpolation.
        Pointers, vectors or files can all be used as input.
        We can also use it to calculate derivatives.
        The interpolation is done using gsl functions.

        \~german \brief Wrapper-Klasse für gsl 1D-Interpolation
    
        Diese Klasse bietet die Möglichkeit eine 1D-Interpolation durchzuführen.
        Dabei können als Eingaben pointer, Vektoren oder ein File verwendet werden.
        Zudem besteht die Möglichkeit sich die Ableitung zurückgeben zu lassen.
        Interpoliert wird mit gsl-Funktionen
    */
	class Interpolation1D{
        	public:

            		//! \~english constructor (does nothing) \~german Konstruktor (macht nix)
            		Interpolation1D();
            
		    /*! \~english \brief constructor (pointer-input)
		        
		        Create a new Interpolation with pointer-input
		        @param x_array           nescessary condition x_array[i+1]>x_array[i]
		        @param y_array           values to be interpolated
		        @param array_length      length of x_array = length of y_array
		        \sa Interpolation1D(std::vector<double> x_array, std::vector<double> y_array)
		        \sa Interpolation1D(std::string filename, uint x_column, uint y_column)

		        \~german \brief Konstruktor (Pointer-Input)
		        
		        Create a new Interpolation with pointer-input
		        @param x_array           Es muss gelten x_array[i+1]>x_array[i]
		        @param y_array           Werte, die interpoliert werden sollen
		        @param array_length      Länge des x_arrays = Länge des y_arrays
		        \sa Interpolation1D(std::vector<double> x_array, std::vector<double> y_array)
		        \sa Interpolation1D(std::string filename, uint x_column, uint y_column)
		    */
		    Interpolation1D(double* x_array, double* y_array, uint array_length);
		    
		    /*! \~english \brief constructor (vector-input)

		        Create a new Interpolation with vector-input
		        @param x_array           nescessary condition x_array[i+1]>x_array[i]
		        @param y_array           values to be interpolated
		        \sa Interpolation1D(double* x_array, double* y_array, uint array_length)
		        \sa Interpolation1D(std::string filename, uint x_column, uint y_column)

		        \~german \brief Konstruktor (Vektor-Input)

		        Create a new Interpolation with vector-input
		        @param x_array           Es muss gelten x_array[i+1]>x_array[i]
		        @param y_array           Werte, die interpoliert werden sollen
		        \sa Interpolation1D(double* x_array, double* y_array, uint array_length)
		        \sa Interpolation1D(std::string filename, uint x_column, uint y_column)
		    */
		    Interpolation1D(std::vector<double> x_array, std::vector<double> y_array);
		    
		    /*! \~english \brief constructor (file-input)
		        Create a new Interpolation with file-input
		        @param filename          file containing nescessary values as columns (separated by tabs)
		        @param x_column          number of column containing x_array
		        @param y_column          number of column containing y_arrays
		        \sa Interpolation1D(double* x_array, double* y_array, uint array_length)
		        \sa Interpolation1D(std::vector<double> x_array, std::vector<double> y_array)

		        \~german \brief Konstruktor (File-Input)

		        Create a new Interpolation with file-input
		        @param filename          Datei, in der die benötigten Werte in Spalten stehen (getrennt durch tabs)
		        @param x_column          Spalte des x_arrays
		        @param y_column          Spalte des y_arrays
		        \sa Interpolation1D(double* x_array, double* y_array, uint array_length)
		        \sa Interpolation1D(std::vector<double> x_array, std::vector<double> y_array)
		    */
		    Interpolation1D(std::string filename, uint x_column, uint y_column);
		    
		    //! \~english destructor \~german Destruktor
		    ~Interpolation1D();
		    
		    /*! \~english \brief returns interpolated value
		        \~german  \brief Rückgabe eines interpolierten Wertes 

		        \~
		        @param xvalue
		        @return f(x) 
		    */
		    double getValue(double xvalue) const;
		    
		    /*! \~english \brief returns derivative at certain value
		        \~german \brief Rückgabe der Ableitung an einer Stelle

		        \~
		        @param xvalue
		        @return f'(x) 
		    */
		    double getDeriv(double xvalue) const;
		    
		    /*! \~english \brief returns 2. derivative at certain value
		        \~german \brief Rückgabe der 2. Ableitung an einer Stelle

		        \~
		        @param xvalue
		        @return f''(x) 
		    */
		    double getDeriv2(double xvalue) const;
		    
		    /*! \~english \brief returns smallest x-value for the interpolation
		        \~german \brief Rückgabe des minimalen x-Wertes für die Interpolation

		        \~
		        @return minimal x value 
		    */
		    double getMin() const;
		    
		    /*! \~english \brief returns largest x-value for the interpolation
		        \~german \brief Rückgabe des maximalen x-Wertes für die Interpolation
		        
		        \~
		        @return maximal x value 
		    */
		    double getMax() const;
		    
		    /*! \~english \brief returns length of input data
		        \~german \brief Rückgabe der Länge der Inputdaten

		        \~
		        @return array_length_ 
		    */
		    uint getLength() const;
		    
		    /*! \~english \brief returns input x-array
		        \~german \brief Rückgabe des Input-x-Arrays

		        \~
		        @return x_array_
		    */
		    double* getXArray() const;
		    
		    /*! \~english \brief returns input y-array
		        \~german \brief Rückgabe des Input-y-Arrays

		        \~
		        @return y_array_ 
		    */
		    double* getYArray() const;
		    
		    /*! \~english \brief Update the y Values of this interpolation with pointer input
		        \~german \brief Update the y Values of this interpolation with pointer input

		        \~
		        @param y_array
		        \sa update(std::vector<double> y_array)
		    */
		    void update(double* y_array);
		    
		    /*! \~english \brief Update the y Values of this interpolation with vector input
		        \~german \brief Update the y Values of this interpolation with vector input

		        \~
		        @param y_array
		        \sa update(double* y_array)
		    */
		    void update(std::vector<double> y_array);
		
		
		private:
		    /*! \~english \brief starts interpolation process
		    
		        initializes interpolation and calls \a update();
		        called by constructors

		        \~german \brief startet den Interpolations-Vorgang
		    
		        initialisiert die Interpolation und ruft \a update() auf;
		        wird von den Konstruktoren aufgerufen

		        \~
		        @param yarray
		    */
		    void createInterpolation(double* yarray);

		    //! \~english saves x-values of the interpolation \~german speichert die x-Werte der Interpolation 
		    double* x_array_;
		    
		    //! \~english saves y-values of the interpolation \~german speichert die y-Werte der Interpolation
		    double* y_array_;

		    //! \~english saves length of input data \~german speichert die Länge der Inputdaten 
		    uint array_length_;
		    
		    //! \~ gsl object
		    gsl_interp_accel* acclerator_;
		    //! \~ gsl object
		    gsl_spline* spline_;
	};
}
#endif /* INTERPOLATION1D_H */
