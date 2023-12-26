#include <iostream>
#include <cmath>
#include <array>
#include <numbers>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "main.h"
#include <pybind11/numpy.h>


namespace py = pybind11;
using namespace std;



double norm(std::array<double, 3Ui64> &r_diff) {
    return std::sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1] + r_diff[2]*r_diff[2]);
}

double norm_np_sq(py::array_t<double> vector_numpy){
    py::buffer_info res_buf = vector_numpy.request();
    double *vector = static_cast<double *>(res_buf.ptr);
    return std::sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
}
double I1(double z, double x_sqr_y_sqr, double r_diff_norm) {
    return z / (x_sqr_y_sqr * r_diff_norm);
}

double I2(double z, double x_sqry_sqr, double r_diff_norm) {

    double num = z * ( 3 * x_sqry_sqr + 2 * z*z);
    double den = 3 * x_sqry_sqr * x_sqry_sqr *r_diff_norm*r_diff_norm*r_diff_norm;
    return num / den;
}

double I3(double z, double x_sqry_sqr, double r_diff_norm, double z2){
    double num = (
            (-2 * (x_sqry_sqr*x_sqry_sqr) * z2)
            + (x_sqry_sqr * (z*z*z + 3 * (z2*z2) * z))
            - 2 * (z2*z2) * -(z*z*z)
        );

    double den = 3 * x_sqry_sqr*x_sqry_sqr * r_diff_norm*r_diff_norm*r_diff_norm;
    return num / den;
}

double lI1(double x_sqry_sqr){
    return 1/x_sqry_sqr;
}

double lI2(double x_sqry_sqr){
    return 2/(3*x_sqry_sqr*x_sqry_sqr);
}

double lI3(double x_sqry_sqr, double z2){
    return (x_sqry_sqr + 2*z2*z2)/(3*x_sqry_sqr*x_sqry_sqr);
}
