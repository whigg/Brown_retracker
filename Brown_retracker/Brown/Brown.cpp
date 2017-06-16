#include "Brown.h"
#include "misc/misc.h"
#include "retrack_func/retracker_functions.h"

// #include <levmar.h>
#include <lib/levmar/levmar.h>
#include <vector>
#include <iostream>
#include <cassert>

using namespace std;

// default constructor
Brown::Brown(int startIndexThermalNoise, int endIndexThermalNoise, double scaling,
                     double timeResolution, vector<double> brownParamInitGuess,
                     double widthRadarPointTargetResponse, bool ThermalNoiseComputIsMedian,
                     double satHeight, double satLat,
                     double antennaBeamWidth) :
    m_c(299792458),
    m_ThermalNoiseIsMedian(ThermalNoiseComputIsMedian)
{
    set_m_r_t(timeResolution);
    set_brown_param_init_guess(brownParamInitGuess);
    set_m_sigma_p(widthRadarPointTargetResponse);
    set_m_scaling(scaling);
    set_antennaBeamWidth(antennaBeamWidth);
    set_startEndIndexThermalNoise(startIndexThermalNoise, endIndexThermalNoise);
    set_satHeight_satLat(satHeight, satLat);
    set_m_a();

}


void Brown::print_config() const
{
    cout << "start index (included) thermal noise computation: " << m_start << " (unit: none)" << endl;
    cout << "end index (excluded) thermal noise computation: " << m_end << " (unit: none)" << endl;
    cout << "time resolution: " << m_r_t << " (unit: s)" << endl;
    cout << "width radar point target response: " << m_sigma_p << " (unit: s)" << endl;
    cout << "thermal noise: " << m_T_n << " (unit: power_echo)" << endl;
    cout << "satellite height: " << m_h << " (unit: m)" << endl;
    cout << "earth radius: " << m_R_e << " (unit: m)" << endl;
    // cout << "skewness: " << m_lambda_s << " (unit: none)" << endl;
    cout << "antenna_beam_width: " << m_theta_0 << " (unit: m)" << endl;
    cout << "gamma: " << m_gamma << " (unit: none)" << endl;
    cout << "a: " << m_a << " (unit: s^-1)" << endl;
    cout << "scaling: " << m_scaling << " (unit: none)" << endl;
    cout << "speed of light: " << m_c << " (unit: m/s)" << endl;

    cout << "initial guesses for Brown parameters: " << endl;
    cout << "   off-nadir mispointing angle: " << m_brown_param_init_guess.at(0) << " (unit: rad)" << endl;
    cout << "   amplitude of the signal: " << m_brown_param_init_guess.at(1) << " (unit: power_echo)" << endl;
    cout << "   epoch: " << m_brown_param_init_guess.at(2) << " (unit: bins)" << endl;
    cout << "   slope of the leading edge: " << m_brown_param_init_guess.at(3) << " (unit: s)" << endl;
}

// =================
// get methods
// =================
int Brown::get_startIndexThermalNoise() const {return m_start;}
int Brown::get_endIndexThermalNoise() const {return m_end;}
int Brown::get_speedOfLight() const {return m_c;}
double Brown::get_scaling() const {return m_scaling;}
double Brown::get_timeResolution() const {return m_r_t;}
vector<double> Brown::get_dataPoints() const {return m_data_points;}
vector<double> Brown::get_brownParamInitGuess() const {return m_brown_param_init_guess;}
double Brown::get_widthRadarPointTargetResponse() const {return m_sigma_p;}
double Brown::get_thermalNoise() const {return m_T_n;}
double Brown::get_satelliteHeight() const {return m_h;}
double Brown::get_EarthRadius() const {return m_R_e;}
double Brown::get_antennaBeamWidth() const {return m_theta_0;}
double Brown::get_gamma() const {return m_gamma;}
double Brown::get_a() const {return m_a;}
vector<double> Brown::get_fit_parameters_initial() const {return m_brown_param_init_guess;}
vector<double> Brown::get_fit_parameters_estimated() const {return m_brown_param_estimated;}

map<string, double> Brown::get_config() const
{

    map<string, double> config;

    // store the config parameters
    config.insert(std::make_pair("start", m_start));
    config.insert(std::make_pair("end", m_end));
    config.insert(std::make_pair("time_resolution", m_r_t));
    config.insert(std::make_pair("width_radar_point_target_response", m_sigma_p));
    config.insert(std::make_pair("thermal_noise", m_T_n_scaled));
    config.insert(std::make_pair("satellite_height", m_h));
    config.insert(std::make_pair("earth_radius", m_R_e));
    config.insert(std::make_pair("skewness", m_lambda_s));
    config.insert(std::make_pair("antenna_beam_width", m_theta_0));
    config.insert(std::make_pair("gamma", m_gamma));
    config.insert(std::make_pair("a", m_a));
    config.insert(std::make_pair("scaling", m_scaling));

    return config;
}


double Brown::fit_value(double t) const
{


    // brown parameters
    vector <double> p;
    p = m_brown_param_estimated;

    // units conversion
    t = t * m_r_t;  // converts the time delay from bins to seconds
    p[2] = p[2] * m_r_t;


    //------ Compute intermediate variables ------//
    double a_ksi = exp((-4*pow(sin(p[0]), 2.0)) / m_gamma);
    double b_ksi = cos(2*p[0]) - ( (pow(sin(2*p[0]), 2.0)) / m_gamma );
    double c_ksi = b_ksi * m_a;

    // Need the current x, i.e. the time from the for loop and/or a parameter that is estimated
    double sigma_c_2 = pow(m_sigma_p, 2.0) + pow(p[3], 2.0);
    double sigma_c = sqrt(sigma_c_2);

    //------ Intermediate equations ------//
    double v = c_ksi * (t - p[2] - 0.5*c_ksi*sigma_c_2);
    double u = (t - p[2] - c_ksi*sigma_c_2) / ( sqrt(2)*sigma_c );

    //------ Forking Brown vs Hayne implementation ------//
    if (m_lambda_s == 0){  // Brown
        return m_T_n + 0.5*a_ksi*p[1]*exp(-v) * (1 + erf(u));

    } else {  // Hayne equation
        double erf_u = erf(u);  // only needed to speed Hayne model (because called twice in the function)
        double part1 = 0.5*a_ksi*p[1]*exp(-v);
        double part2 = (1+erf_u);
        double part3 = (m_lambda_s/6) * (pow(p[3], 2.0)/sigma_c_2);
        double part4 = part2 * pow(c_ksi, 3.0)*pow(sigma_c, 3.0);
        double part5 = (sqrt(2)/sqrt(M_PI));
        double part6 = 2*pow(u, 2.0);
        double part7 = 3*sqrt(2)*c_ksi*sigma_c*u;
        double part8 = 3*pow(c_ksi ,2.0)*pow(sigma_c, 2.0) - 1;
        double part9 = exp(-pow(u, 2.0));

        // double test_point =  m_T_n + part1 * (part2 + part3 * (part4 - (part5 * (part6 + part7 + part8) * part9)));

        return m_T_n + part1 * (part2 + part3 * (part4 - (part5 * (part6 + part7 + part8) * part9)));
    }

}


vector<double> Brown::get_fit_values() const
{

    unsigned int dataPointCount = m_data_points.size();
    vector<double> fitted_values(dataPointCount, -9999);
    for (unsigned int bin=0; bin < dataPointCount; bin++) fitted_values[bin] = fit_value(bin);

    return fitted_values;

}

// =================
// set methods
// =================
void Brown::set_startEndIndexThermalNoise(int startIndexThermalNoise, int endIndexThermalNoise)
{

    // check inputs
//    assert(startIndexThermalNoise > 0);
//    assert(startIndexThermalNoise < endIndexThermalNoise);
    if (startIndexThermalNoise < 0) throw std::range_error( "'startIndexThermalNoise' cannot be negative\n" );
    if (startIndexThermalNoise >= endIndexThermalNoise ) throw range_error( "'startIndexThermalNoise' must be smaller than 'endIndexThermalNoise'\n" );
    if ((endIndexThermalNoise - startIndexThermalNoise) < 1) throw range_error( "end-start (i.e. subset of the vector to compute the thermal noise with) must have at least 1 element \n" );

    m_start = startIndexThermalNoise;
    m_end = endIndexThermalNoise;

}


void Brown::set_thermalNoiseIsMedian(bool is_median)
{
    m_ThermalNoiseIsMedian = is_median;

    // if already datavector, compute the new thermal noise
    // if (){}
}



void Brown::set_thermalNoise()
{

    //  check start and end index are within the vector elements
    if (m_end > m_data_points.size()) throw invalid_argument( "'m_end' must be smaller than 'm_data_points' length" );

    // compute the thermal noise with mean or median
    if (m_ThermalNoiseIsMedian){
        m_T_n = median(m_data_points, m_start, m_end);
    } else {
        m_T_n = mean(m_data_points, m_start, m_end);
    }


}



//void Brown::set_dataPoints(vector<double> dataPoints)
//{
//    // store the waveform datapoints
//    m_data_points = dataPoints;

//    // compute the thermal noise
//    set_thermalNoise();

//}

void Brown::set_m_a()
{
    m_a = 4*m_c / ( m_gamma * m_h * (1 + (m_h / m_R_e)) );
}



void Brown::set_m_R_e()
{

    double lat_rad = m_lat;

    // compute the earth radius
    const double r1 = 6378137;  // equatorial earth radius
    const double r2 = 6356752.3142;  // polar earth radius

    double part1 = pow(cos(lat_rad), 2.0) / pow(r1, 2.0);
    double part2 = pow(sin(lat_rad), 2.0) / pow(r2, 2.0);
    double R_e = pow(part1 + part2, -0.5);  // radius of the earth at the specific lat

    // set the earth radius
    m_R_e = R_e;

}

void Brown::set_satHeight_satLat(double satHeight, double satLat)
{

    // check satHeight is plausible
    if (satHeight <= 0) throw invalid_argument("'satelliteHeight' must be > 0.\n");

    // check satLat is plausible
    double lat_deg = rad2deg(satLat);  // converts in degrees
    if ((lat_deg < -180) || (lat_deg > 180)) {
        throw invalid_argument( "'satLat' (i.e. latitude in deg) must be in [-180;180].\n" );
    }

    // set satHeight, satLat and unit
    m_h = satHeight;
    m_lat = satLat;

    // estimate and set earth radius
    set_m_R_e();

}


void Brown::set_m_scaling(double scaling)
{
    // check scaling makes sens (i.e. positive)
    if (scaling <= 0) throw range_error("'scaling' must be strictely positive. \n");
    m_scaling = scaling;
}

//double Brown::set_widthRadarPointTargetResponse() {return m_sigma_p;}
//double Brown::set_thermalNoise() {return m_T_n;}
//double Brown::set_satelliteHeight() {return m_h;}
//double Brown::set_EarthRadius() {return m_R_e;}
//double Brown::set_skewness() {return m_lambda_s;}
//double Brown::set_gamma() {return m_gamma;}
//double Brown::set_a() {return m_a;}


void Brown::set_m_r_t(double timeResolution)
{
    // check input variable
    if (timeResolution <= 0) throw range_error("'timeResolution' must be strictely positive. \n");

    // set the value
    m_r_t = timeResolution / m_c;

}

void Brown::set_m_sigma_p(double widthRadarPointTargetResponse)
{
    // check valid input
    if (widthRadarPointTargetResponse <= 0) throw range_error("'widthRadarPointTargetResponse' must be strictely positive. \n");

    // convert to seconds
    m_sigma_p = widthRadarPointTargetResponse * m_r_t;
}

void Brown::set_antennaBeamWidth(double antennaBeamWidth)
{
    // check input is valid
    if (antennaBeamWidth <= 0) throw range_error("'antennaBeamWidth' must be strictely positive. \n");

    // set the value
    m_theta_0 = antennaBeamWidth;

    // compute gamma
    m_gamma = pow(sin(m_theta_0), 2.0) / (2*log(2));

}

void Brown::set_brown_param_init_guess(std::vector<double> brownParamInitGuess)
{
    // check input parameters are valid
    if (brownParamInitGuess.size() != 4) throw range_error("'brownParamInitGuess' must contain 4 elements. \n");
    if (brownParamInitGuess.at(0) < 0) throw range_error("'brownParamInitGuess.at(0)', i.e. the off-nadir mispointing angle, must be positive. \n");
    if (brownParamInitGuess.at(1) < 0) throw range_error("'brownParamInitGuess.at(1)', i.e. the amplitude of the signal (in unit of the waveform power echoe)"
                                                ", must be positive. \n");
    if (brownParamInitGuess.at(2) < 0) throw range_error("'brownParamInitGuess.at(2)', i.e. the epoch (in bins), must be positive. \n");
    if (brownParamInitGuess.at(3) < 0) throw range_error("'brownParamInitGuess.at(3)', i.e. the slope of the leading edge (in seconds), must be positive. \n");

    // store the model parameters initial guesses
    m_brown_param_init_guess = brownParamInitGuess;

}

void Brown::fit_brown(std::vector<double> data, double satHeight, double satLat)
{
    fit(data, satHeight, satLat, 0);
}

void Brown::fit_hayne(std::vector<double> data, double satHeight, double satLat, double skewness)
{
    fit(data, satHeight, satLat, skewness);
}

void Brown::fit(std::vector<double> data, double satHeight, double satLat, double skewness)
{

    // Set-up
    // ======

    m_lambda_s = skewness;

    // Max number of iterations in the optimization problem
    unsigned int max_iter = 1000;

    // store data
    m_data_points = data;

    // Dimensionality of the problem
    const unsigned int brownParameterCount = 4;
    unsigned int dataPointCount = m_data_points.size();

    // check and store satHeight and lat and compute earth radius a lat
    set_satHeight_satLat(satHeight, satLat);

    // scale the data points
    m_data_points_scaled = m_data_points;
    double data_points_scaled[dataPointCount];
    for (unsigned int pnt = 0; pnt < m_data_points_scaled.size(); pnt++) {
        m_data_points_scaled[pnt] = m_data_points_scaled[pnt] * m_scaling;
        data_points_scaled[pnt] = m_data_points_scaled[pnt];
    }

    // Compute the thermal noise
    set_thermalNoise();  // unscaled
    m_T_n_scaled = m_T_n * m_scaling;  // scaled

    // Create the map containing the info needed by the levmar as an external fit
    std::map<std::string, double> config;
    config = get_config();


    // make a pointer out of the map params to make it possible to pass them to levmar algo
    void *config_void_ptr = &config;

    // converting param vect to array with changed units for the levmar algo
    double brown_param_init_guess_converted[brownParameterCount];
    for (unsigned int param=0; param<brownParameterCount; param++) {
        brown_param_init_guess_converted[param] = m_brown_param_init_guess[param];
    }

    brown_param_init_guess_converted[1] = brown_param_init_guess_converted[1] * m_scaling;
    brown_param_init_guess_converted[2] = brown_param_init_guess_converted[2] * m_r_t;

    /* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
    double opts[LM_OPTS_SZ];  // create a 4 or 5 elements array for levmar options
    opts[0]=LM_INIT_MU;     // tau: scalefactor for initial mu
    opts[1]=1E-15;          // epsilon1: stopping threshold for ||J^T e||_inf
    opts[2]=1E-32;          // epsilon2: stopping threshold for ||Dp||_2
    opts[3]=1E-20;          // epsilon3: stopping threshold for ||e||_2
    opts[4]=1E-30;  // delta: step used in approximation of Jacobian. Relevant only if the finite difference Jacobian version is used
//    opts[1]=1E-60;          // epsilon1: stopping threshold for ||J^T e||_inf
//    opts[2]=1E-60;          // epsilon2: stopping threshold for ||Dp||_2
//    opts[3]=1E-60;          // epsilon3: stopping threshold for ||e||_2
//    opts[4]=1E-60;  // delta: step used in approximation of Jacobian. Relevant only if the finite difference Jacobian version is used


    double brownInfo[LM_INFO_SZ];

    // Fit
    // ===
    dlevmar_dif(brown, brown_param_init_guess_converted, data_points_scaled,
                              brownParameterCount, dataPointCount, max_iter, opts,
                              brownInfo, NULL, NULL, config_void_ptr);


    // Store the parameters
    // ====================

    // info regarding the fit procedure
    const char *brown_info_names[] = {"||e||_2_init", "||e||_2", "||J^T e||_inf",
                                      "||Dp||_2", "mu_max[J^T J]_ii", "iterations",
                                      "terminating_reason", "func_evals", "jacobian_evals", "lin_systems_solved"};
    for (unsigned int i=0; i<(sizeof(brownInfo)/sizeof(*brownInfo)); i++)
    {
        m_info_fit.insert(std::make_pair(brown_info_names[i], brownInfo[i]));
    }

    // estimated model parameters rescaled
    vector<double> brown_param_estimated(brownParameterCount, -9999);
    for (unsigned int param=0; param<brownParameterCount; param++) {
        brown_param_estimated[param] = brown_param_init_guess_converted[param];
    }

    m_brown_param_estimated = brown_param_estimated;
    m_brown_param_estimated[1] = m_brown_param_estimated[1] / m_scaling;
    m_brown_param_estimated[2] = m_brown_param_estimated[2] / m_r_t;

    // info regarding the options of the fit procedure (perhaps allow user to modify them later)
    const char *opt_names[] = {"scale_factor_init_mu", "stop_threshold_for_||J^T e||_inf",
                                      "stop_threshold_for_||Dp||_2", "stop_threshold_for_||e||_2",
                                      "step_size_for_difference_approx_of_Jacobian"};
    for (unsigned int i=0; i<(sizeof(opts)/sizeof(*opts)); i++){
        m_opts.insert(std::make_pair(opt_names[i], opts[i]));
    }


}
