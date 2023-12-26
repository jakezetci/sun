#include <iostream>
#include <cmath>
#include <array>
#include <numbers>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "main.h"
#include <pybind11/numpy.h>
#include <time.h>
#include <ctime>
#include <chrono>
#include "lib.h"

class Timer
{
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const { 
        return std::chrono::duration_cast<second_>
            (clock_::now() - beg_).count(); }

private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace std;

int add(int i, int j) {
    return i + j;
}


std::array<double,3> GreenBl(double r1x, double r1y, double r1z, double r2x, double r2y, double r2z)
{
    std::array<double,3> G, r1, r2, r_diff;
    r1[0] = r1x;
    r1[1] = r1y;
    r1[2] = r1z;
    r2[0] = r2x;
    r2[1] = r2y;
    r2[2] = r2z;
    double a = 696340000;
    
    for(int i = 0 ; i < 3; i++){
        r_diff[i] = r1[i] - r2[i];
    }

    /* preparing things for I[i] functions computation
    so we don't do same things multiple times*/
    double r_diff_norm = norm(r_diff);
    double x_sqr_y_sqr = r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1];
    double x = r_diff[0];
    double y = r_diff[1];
    double z = r_diff[2];
    // computing the functions
    double ili1, ili2, ili3;    
    
    ili1 = I1(z, x_sqr_y_sqr, r_diff_norm) - lI1(x_sqr_y_sqr);
    ili2 = I2(z, x_sqr_y_sqr, r_diff_norm) - lI2(x_sqr_y_sqr);
    ili3 = I3(z, x_sqr_y_sqr, r_diff_norm, r2[2]) - lI3(x_sqr_y_sqr, r2[2]);

    double xya = (r1[0]*r1[0] + r1[1]*r1[1] - a*a)*ili2 +ili3;
    double den = 4 *std::numbers::pi * a;
    G[0] = (2 * r1[0] * ili1 - 3 * x * xya)/den;
    G[1] = (2 * r1[1] * ili1 - 3 * y * xya)/den;
    G[2] = ((r1[0]*r1[0] + r1[1]*r1[1]+ r1[2]*r1[2] - a*a) / (r_diff_norm*r_diff_norm*r_diff_norm))/den;
    
    return G;
}




py::array_t<double> B_comp(py::array_t<double> r, py::array_t<double> values, 
    py::array_t<double> points,py::array_t<double> areas) {
    py::buffer_info r_buf = r.request(), pts_buf = points.request(), values_buf = values.request(), areas_buf = areas.request();
    py::array_t<double> result = py::array_t<double>({3});
    py::buffer_info res_buf = result.request();
    double *ptr_res = static_cast<double *>(res_buf.ptr);
    ptr_res[0] = 0.0;
    ptr_res[1] = 0.0;
    ptr_res[2] = 0.0;
    
    int N = static_cast<int>(values_buf.shape[0]);
    double *ptr_pts = static_cast<double *>(pts_buf.ptr);
    double *ptr_values =  static_cast<double *>(values_buf.ptr);
    double *ptr_areas = static_cast<double *>(areas_buf.ptr);
    double *ptr_r = static_cast<double *>(r_buf.ptr);
    double r_x = ptr_r[0];
    double r_y = ptr_r[1];
    double r_z = ptr_r[2];
    double multiplier;
    std::array< double,3 > G;
    for(int i = 0; i < N; i++){
        multiplier = ptr_values[i]*ptr_areas[i];
        G = GreenBl(r_x, r_y, r_z, ptr_pts[i*3], ptr_pts[i*3+1], ptr_pts[i*3+2]);
        for (int j = 0; j < 3; j++){
            ptr_res[j] +=  multiplier * G[j];
        };
    }

    return result;
}



double bitmap_energy(py::array_t<double> xyz, double volume, py::array_t<double> values, 
    py::array_t<double> points,py::array_t<double> areas, int timestamp = 10){
    py::buffer_info xyz_buf = xyz.request();
    double *ptr_xyz = static_cast<double *>(xyz_buf.ptr);
    int N = static_cast<int>(xyz_buf.shape[0]);

    py::array_t<double> r = py::array_t<double>({3});
    py::buffer_info r_buf = r.request();
    double *ptr_r = static_cast<double *>(r_buf.ptr);
    double energy = 0.0;
    py::array_t<double> B = py::array_t<double>({3});
    py::buffer_info B_buf = B.request();

    Timer tmr;
    
    for(int i = 0; i < N; i++){
        ptr_r[0] = ptr_xyz[i*3];
        ptr_r[1] = ptr_xyz[i*3+1];
        ptr_r[2] = ptr_xyz[i*3+2];
        B = B_comp(r, values, points, areas);
        double B_norm = norm_np_sq(B);
        energy += B_norm * volume /(8*std::numbers::pi);
        if (i % timestamp == 0){
            if (i != 0){
                double t = tmr.elapsed();
                std::cout << "Values " << i-timestamp << "-"<< i << " done in ";
                std::cout << t << " seconds" << std::endl;
                tmr.reset();
            };
        
        };
    };
    return energy;
}

void print(std::array<double, 3Ui64> &r1)
{
    for (int i = 0; i < 3; i++)
    {
        std::cout << r1[i] << " ";
    }
}
int multiply(int i, int j) {
    return i * j;
}


PYBIND11_MODULE(cpp_module, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: cpp_module

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.def("multiply", &multiply, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.def("green", &GreenBl, R"pbdoc(
        Computes Green's function
    )pbdoc");

    

    m.def("b_comp", &B_comp, "computes the magnetic field");
    m.def("energy", &bitmap_energy, "gaga");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
