#include "retracker_functions.h"
#include <cmath>
#include <math.h>
#include <iostream>
#include <stdexcept>  // needed to throw exceptions
#include "misc/misc.h"
#include <map>
#include <string>
#include <iterator>
using namespace std;




/*==============================================================================*/
/* Full Brown model */


//------------------ to optimize parameters -------------------------//
void brown(double *p, double *y, int , int n, void *config_void_ptr)
{

    // cast the void pointer to a map pointer (needed before being able to dereference it)
    std::map<std::string, double> *config_map_ptr = static_cast<std::map<std::string, double>*>(config_void_ptr);
    // dereferencing the map pointer to get back the original map containing the config parameters
    std::map<std::string, double> config = *config_map_ptr;

    // Extract the config parameters from the map container
    double r_t = config.at("time_resolution");
    double sigma_p = config.at("width_radar_point_target_response");
    double T_n = config.at("thermal_noise");
    double lambda_s = config.at("skewness");
    double gamma = config.at("gamma");
    double a = config.at("a");
    double scaling = config.at("scaling");

    // initialization of some variables
    double t;  // time delay
    double a_ksi;
    double b_ksi;
    double c_ksi;
    double sigma_c_2;
    double sigma_c;
    double v;
    double u;
    double erf_u;
    double part1;
    double part2;
    double part3;
    double part4;
    double part5;
    double part6;
    double part7;
    double part8;
    double part9;






    //------ Compute the brown function ------//

    for(int i=0; i<n; i++){

        //------ Compute intermediate variables ------//

        a_ksi = exp((-4*pow(sin(p[0]), 2.0)) / gamma);
        b_ksi = cos(2*p[0]) - ( (pow(sin(2*p[0]), 2.0)) / gamma );
        c_ksi = b_ksi * a;

        // Need the current x, i.e. the time from the for loop and/or a parameter that is estimated
        sigma_c_2 = pow(sigma_p, 2.0) + pow(p[3], 2.0);
        sigma_c = sqrt(sigma_c_2);


        t = i * r_t;  // get time in seconds from bin nb i and time resolution
        v  = c_ksi * (t - p[2] - (c_ksi*sigma_c_2 / 2));
        u = (t - p[2] - c_ksi*sigma_c_2) / ( sqrt(2)*sigma_c );


        if (lambda_s == 0){  //------ Final Brown equation ------//
            y[i] = T_n + 0.5*a_ksi*p[1]*exp(-v) * (1 + erf(u));

        } else {  //------ Final Hayne equation ------//

            // Hayne model
            erf_u = erf(u);  // only needed to speed Hayne model (because called twice in the function)
            part1 = 0.5*a_ksi*p[1]*exp(-v);
            part2 = (1+erf_u);
            part3 = (lambda_s/6) * (pow(p[3], 2.0)/sigma_c_2);
            part4 = part2 * pow(c_ksi, 3.0)*pow(sigma_c, 3.0);
            part5 = (sqrt(2)/sqrt(M_PI));
            part6 = 2*pow(u, 2.0);
            part7 = 3*sqrt(2)*c_ksi*sigma_c*u;
            part8 = 3*pow(c_ksi ,2.0)*pow(sigma_c, 2.0) - 1;
            part9 = exp(-pow(u, 2.0));

            y[i] = T_n + part1 * (part2 + part3 * (part4 - (part5 * (part6 + part7 + part8) * part9)));
        }

    }

}


//------------------ to plot -------------------------//
double brown(double *p, double t, std::map<std::string, double> config)
{

    // extract the config parameters
    double T_n = config.at("thermal_noise");  // unit: power unit
    double r_t = config.at("time_resolution");  // unit: s
    double sigma_p = config.at("width_radar_point_target_response");  // unit: s
    double gamma = config.at("gamma");
    double a = config.at("a");
    double scaling = config.at("scaling");
    double lambda_s = config.at("skewness");  // unit: none

    // descale
    T_n = T_n / scaling;
    double p1 = p[1] / scaling;


    // units conversion
    t = t * r_t;  // converts the time delay from bins to seconds

    //------ Compute intermediate variables ------//
    double a_ksi = exp((-4*pow(sin(p[0]), 2.0)) / gamma);
    double b_ksi = cos(2*p[0]) - ( (pow(sin(2*p[0]), 2.0)) / gamma );
    double c_ksi = b_ksi * a;

    // Need the current x, i.e. the time from the for loop and/or a parameter that is estimated
    double sigma_c_2 = pow(sigma_p, 2.0) + pow(p[3], 2.0);
    double sigma_c = sqrt(sigma_c_2);

    //------ Intermediate equations ------//
    double v  = c_ksi * (t - p[2] - 0.5*c_ksi*sigma_c_2);
    double u = (t - p[2] - c_ksi*sigma_c_2) / ( sqrt(2)*sigma_c );

    //------ Forking Brown vs Hayne implementation ------//
    if (lambda_s == 0){  // Brown
        return T_n + 0.5*a_ksi*p1*exp(-v) * (1 + erf(u));

    } else {  // Hayne equation
        double erf_u = erf(u);  // only needed to speed Hayne model (because called twice in the function)
        double part1 = 0.5*a_ksi*p1*exp(-v);
        double part2 = (1+erf_u);
        double part3 = (lambda_s/6) * (pow(p[3], 2.0)/sigma_c_2);
        double part4 = part2 * pow(c_ksi, 3.0)*pow(sigma_c, 3.0);
        double part5 = (sqrt(2)/sqrt(M_PI));
        double part6 = 2*pow(u, 2.0);
        double part7 = 3*sqrt(2)*c_ksi*sigma_c*u;
        double part8 = 3*pow(c_ksi ,2.0)*pow(sigma_c, 2.0) - 1;
        double part9 = exp(-pow(u, 2.0));

        // double test_point =  T_n + part1 * (part2 + part3 * (part4 - (part5 * (part6 + part7 + part8) * part9)));

        return T_n + part1 * (part2 + part3 * (part4 - (part5 * (part6 + part7 + part8) * part9)));
    }

}

