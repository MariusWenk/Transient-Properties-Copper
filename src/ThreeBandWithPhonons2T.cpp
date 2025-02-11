/* ************************************
 * File:   ThreeBandWithPhonons2T.hpp
 * Author: ndione
 *
 * Created on June 17, 2020, 5:00 PM
 *************************************/

#include "ThreeBandWithPhonons2T.hpp"

ThreeBandWithPhonons2T::ThreeBandWithPhonons2T() : ODE(9) {

	std::cerr<<"Default constructor of class \"ThreeBandWithPhonons2T\". Please use other constructors\n";
	std::exit(1);
}

ThreeBandWithPhonons2T::ThreeBandWithPhonons2T(std::string sp_filename2D, std::string d_filename2D, Cfermi_distribution *fermi_ptr,  uint col0, uint col1, uint col2) : ODE(16) {
	
}

ThreeBandWithPhonons2T::~ThreeBandWithPhonons2T() {
	
}
