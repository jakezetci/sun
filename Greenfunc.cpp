#include <iostream>
#include <cmath>
#include <array>
using namespace std;

double norm(double r_diff[]) {
    return std::sqrt(r_diff[0]*r_diff[0] + r_diff[1]*r_diff[1] + r_diff[2]*r_diff[2]);
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
    return 1/(3*x_sqry_sqr*x_sqry_sqr);
}

double lI3(double x_sqry_sqr, double z2){
    return (x_sqry_sqr + 2*z2)/(3*x_sqry_sqr*x_sqry_sqr);
}


std::array<double,3> GreenBl(double r1[3], double r2[3])
{
    std::array<double,3> G;
    double a = 696340000;
    double r_diff[3];
    for(int i =1; i < 3; i++){
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
    G[0] = 2 * r1[0] * ili1 - 3 * x * xya;
    G[1] = 2 * r1[1] * ili1 - 3 * y * xya;
    G[2] = (r1[0]*r1[0] + r1[1]*r1[1]+ r1[2]*r1[2] - a*a) / (r_diff_norm*r_diff_norm*r_diff_norm);

    return G;
} 