#include "interpolation1d.hpp"

using namespace ufd;

Interpolation1D::Interpolation1D() {

}

Interpolation1D::Interpolation1D(double* xarray, double* yarray, uint array_length) {
    this->array_length_ = array_length;
    this->x_array_= xarray;
    createInterpolation(yarray);
}

Interpolation1D::Interpolation1D(std::vector<double> xarray, std::vector<double> yarray) {
    this->array_length_ = xarray.size();
    this->x_array_=new double[array_length_];
    for(uint i=0; i<array_length_; i++){
        x_array_[i]=xarray[i];
    }
    createInterpolation(&yarray[0]);
}

Interpolation1D::Interpolation1D(std::string filename, uint x_column, uint y_column) {
    std::vector<double> x_in = BasicInputparser::readColumn(filename, x_column);
    std::vector<double> y_in = BasicInputparser::readColumn(filename, y_column);

    if (x_in.size() != y_in.size()) {
        arError() << "Spalten sind nicht gleich lang";
        exit(1);
    }

    this->array_length_ = x_in.size();
    this->x_array_=new double[array_length_];
    for(uint i=0; i<array_length_; i++){
        x_array_[i]=x_in[i];
    }
    double * yarray = &y_in[0];

    createInterpolation(yarray);
}

Interpolation1D::~Interpolation1D() {
    gsl_spline_free(spline_);
    gsl_interp_accel_free(acclerator_);
}

void Interpolation1D::createInterpolation(double* yarray) {
    acclerator_ = gsl_interp_accel_alloc();
    spline_ = gsl_spline_alloc(gsl_interp_akima, array_length_);
    update(yarray);
}

void Interpolation1D::update(double* y_array) {
    this->y_array_ = y_array;
    gsl_spline_init(spline_, x_array_, y_array_, array_length_);
}

void Interpolation1D::update(std::vector<double> y_array) {
    update(&y_array[0]);
}

double Interpolation1D::getValue(double xvalue) const{
    return gsl_spline_eval(spline_, xvalue, acclerator_);
}

double Interpolation1D::getDeriv(double xvalue) const{
    return gsl_spline_eval_deriv(spline_, xvalue, acclerator_);
}

double Interpolation1D::getDeriv2(double xvalue) const{
    return gsl_spline_eval_deriv2(spline_, xvalue, acclerator_);
}

double Interpolation1D::getMin() const{
    return x_array_[0];
}

double Interpolation1D::getMax() const{
    return x_array_[array_length_-1];
}

uint Interpolation1D::getLength() const{
    return array_length_;
}

double* Interpolation1D::getXArray() const{
    return x_array_;
}

double* Interpolation1D::getYArray() const{
    return y_array_;
}
