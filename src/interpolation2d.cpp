#include "interpolation2d.hpp"

using namespace ufd;

Interpolation2D::Interpolation2D() {

}

Interpolation2D::Interpolation2D(double* x_array, double* y_array, double* z_array, uint x_array_length, uint y_array_length) {
    	this->x_array_length_ = x_array_length;
    	this->y_array_length_ = y_array_length;
    	this->x_array_= x_array;
    	this->y_array_= y_array;
    	createInterpolation(z_array);
}

Interpolation2D::Interpolation2D(double* x_array, double* y_array, double** z_array, uint x_array_length, uint y_array_length) {
    	this->x_array_length_ = x_array_length;
    	this->y_array_length_ = y_array_length;
    	this->x_array_= x_array;
    	this->y_array_= y_array;

    	double ar[x_array_length_*y_array_length_];
    	for (uint i = 0; i < x_array_length_; i++){
        	for (uint j = 0; j < y_array_length_; j++){
            		ar[j*x_array_length_ + i] = z_array[i][j];
        	}
   	}

    	createInterpolation(ar);
}


Interpolation2D::Interpolation2D(std::vector<double> x_array, std::vector<double> y_array, std::vector<double> z_array) {
    	this->x_array_length_ = x_array.size();
    	this->y_array_length_ = y_array.size();
    	this->x_array_=new double[x_array_length_];
    	//\todo Due to problems with references (which I could not yet resolve) the pointer is clumsily copied here
    	for (uint i=0; i<x_array_length_; i++){
        	x_array_[i]=x_array[i];
    	}
    	this->y_array_=new double[y_array_length_];
    	for (uint i=0; i<y_array_length_; i++){
        	y_array_[i]=y_array[i];
    	}

    	createInterpolation(&z_array[0]);
}

Interpolation2D::Interpolation2D(std::vector<double> x_array, std::vector<double> y_array, std::vector<std::vector<double>> z_array) {
    	this->x_array_length_ = x_array.size();
    	this->y_array_length_ = y_array.size();
    	this->x_array_=new double[x_array_length_];
    	//\todo Due to problems with references (which I could not yet resolve) the pointer is clumsily copied here
    	for (uint i=0; i<x_array_length_; i++){
        	x_array_[i]=x_array[i];
    	}
    	this->y_array_=new double[y_array_length_];
    	for (uint i=0; i<y_array_length_; i++){
        	y_array_[i]=y_array[i];
    	}

    	std::vector<double> z_ar(x_array_length_*y_array_length_);
    	for (uint i = 0; i < x_array_length_; i++){
        	for (uint j = 0; j < y_array_length_; j++){
            	z_ar[j*x_array_length_ + i] = z_array[i][j];
        	}
    	}

    	createInterpolation(&z_ar[0]);
}

Interpolation2D::Interpolation2D(std::string filename, uint x_col, uint y_col, uint z_col) {
    	std::vector<double> x_in = ufd::BasicInputparser::readColumn(filename, x_col);
    	std::vector<double> y_in = ufd::BasicInputparser::readColumn(filename, y_col);
    	std::vector<double> z_in = ufd::BasicInputparser::readColumn(filename, z_col);

    	std::vector<double> x_vec;
    	std::vector<double> y_vec;
    	std::vector<double>::iterator it;

    	//find all different x values and collect them in the vector
    	for (uint i = 0; i < x_in.size(); i++){
        	it = std::find(x_vec.begin(),x_vec.end(), x_in[i]);
        	if (it == x_vec.end()){
            		//if value not already in vector add it
            		x_vec.push_back(x_in[i]);
        	}
    	}

    	//find all different y values and collect them in the vector
    	for (uint i = 0; i < y_in.size(); i++){
        	it = std::find(y_vec.begin(),y_vec.end(), y_in[i]);
        	if (it == y_vec.end()){
           	 	//if value not already in vector add it
            		y_vec.push_back(y_in[i]);
        	}
    	}

    	//sort values in vectors ascending
    	std::sort(x_vec.begin(), x_vec.end());
    	std::sort(y_vec.begin(), y_vec.end());

    	/*if (x_vec.size()*y_vec.size() != z_in.size()){
        	arError() << "Length of x,y,z columns does not fit. Exiting...";
        	exit(1);
    	}*/ 

    	std::vector<double> z_vec(x_vec.size() * y_vec.size());
    	std::vector<double>::iterator itx, ity;
    	uint x_idx, y_idx;

    	//put the z values in the correct order in the z_vec
    	for (uint i = 0; i < z_in.size(); i++){
        	itx = std::find(x_vec.begin(), x_vec.end(), x_in[i]);
        	x_idx = distance(x_vec.begin(), itx);
        	ity = std::find(y_vec.begin(), y_vec.end(), y_in[i]);
        	y_idx = distance(y_vec.begin(), ity);
        	z_vec[y_idx * x_vec.size() + x_idx] = z_in[i];
    	}

    	double * zarray = &z_vec[0];
	
    	this->x_array_length_ = x_vec.size();
    	this->y_array_length_ = y_vec.size();
    	this->x_array_=new double[x_array_length_];
    	//Due to problems with references (which I could not yet resolve) the pointer is clumsily copied here
    	for (uint i=0; i<x_array_length_; i++){
        	x_array_[i]=x_vec[i];
    	}
    	this->y_array_=new double[y_array_length_];
    	for (uint i=0; i<y_array_length_; i++){
        	y_array_[i]=y_vec[i];
    	}

    	createInterpolation(zarray);
	delete [] x_array_;
	delete [] y_array_;
}

Interpolation2D::~Interpolation2D() {
    	gsl_spline2d_free(spline_);
    	gsl_interp_accel_free(xacc_);
    	gsl_interp_accel_free(yacc_);
}

void Interpolation2D::createInterpolation(double* z_array) {
    	xacc_ = gsl_interp_accel_alloc();
    	yacc_ = gsl_interp_accel_alloc();
    	spline_ = gsl_spline2d_alloc(T, x_array_length_, y_array_length_);
    	update(z_array);
}

void Interpolation2D::update(std::vector<double> z_array) {
    	update(&z_array[0]);
}

void Interpolation2D::update(std::vector<std::vector<double>> z_array) {
    	std::vector<double> vec(x_array_length_*y_array_length_);

    	for (uint i = 0; i < x_array_length_; i++){
        	for (uint j = 0; j < y_array_length_; j++){
            		vec[j*x_array_length_ + i] = z_array[i][j];
        	}
    	}
    	update(vec);
}

void Interpolation2D::update(double* z_array) {
    	this->z_array_ = z_array;
    	gsl_spline2d_init(spline_, x_array_, y_array_, z_array, x_array_length_, y_array_length_);
}

void Interpolation2D::update(double** z_array){
    	double* ar = new double[x_array_length_*y_array_length_];

    	for (uint i = 0; i < x_array_length_; i++){
        	for (uint j = 0; j < y_array_length_; j++){
            		ar[j*x_array_length_ + i] = z_array[i][j];
        	}
    	}

    	update(ar);
}

double Interpolation2D::getValue(double x_value, double y_value) const{
    	return gsl_spline2d_eval(spline_, x_value, y_value, xacc_, yacc_);
}

/*
   double Interpolation2D::get_deriv(double xvalue, double yvalue) const{
   return gsl_spline2d_eval_deriv(spline_, xvalue, yvalue, xacc_, yacc_);
   }

   double Interpolation2D::get_deriv2(double xvalue, double yvalue) const{
   return gsl_spline2d_eval_deriv2(spline_, xvalue, yvalue, xacc_, yacc_);
   }*/

double Interpolation2D::getMinX() const{
    	return x_array_[0];
}

double Interpolation2D::getMaxX() const{
   	 return x_array_[x_array_length_-1];
}

double Interpolation2D::getMinY() const{
    	return y_array_[0];
}

double Interpolation2D::getMaxY() const{
    	return y_array_[y_array_length_-1];
}

uint Interpolation2D::getXLength() const{
    	return x_array_length_;
}

uint Interpolation2D::getYLength() const{
    	return y_array_length_;
}

double * Interpolation2D::getXArray() const{
    	return x_array_;
}

double * Interpolation2D::getYArray() const{
    	return y_array_;
}

double * Interpolation2D::getZArray() const{
    	return z_array_;
}

