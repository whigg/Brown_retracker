
#include "retrack_func/retracker_functions.h"
#include "misc/misc.h"
#include "Brown/Brown.h"

#include <iterator>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;




//==========================================================================
// Main program

// =============================
// Create a minimal examples
// =============================
int main()
{

    // ====================================
    // Brown fit program on minimal example
    // ====================================

    // load data
    // =========

    // declare variables
    string path_waveforms = "../Brown_retracker_master/data/waveforms_LRM_1b.txt";
    string path_alt = "../Brown_retracker_master/data/alt_LRM_1b.txt";
    string path_lat = "../Brown_retracker_master/data/lat_LRM_1b.txt";
    char delimiter = ';';
    vector<vector<double>> waveforms_full, alts, lats;


    // load the data
    waveforms_full = load_txt(path_waveforms, delimiter);
    alts = load_txt(path_alt, delimiter);
    lats = load_txt(path_lat, delimiter);

    unsigned int N = waveforms_full.size();
    unsigned int D = waveforms_full[0].size();


    // cut the waveform edges (to avoid dealing with the echoes from the previous waveform)
    // ======================

    unsigned int skipBegin=10, skipEnd=10;
    D = D - skipBegin - skipEnd;
    vector<vector<double>> waveforms (N, vector<double> (D, -9999));

    for (unsigned int row=0; row<N; row++) {
        for (unsigned int col=0; col<D; col++) {
            waveforms[row][col] = waveforms_full[row][col+skipBegin];
        }
    }

    // save the cut data
    string path_waveforms_cut = "../Brown_retracker_master/data/outputs/waveforms_LRM_1b_cut.txt";
    save_txt(waveforms, path_waveforms_cut);


    // Create object
    // =============

    // Fast way
    Brown satellite_empty;

    // // Full options setting)
    // Brown satellite_full();


    // Change some parameters at will
    // ==============================

    // Change start end indices
    satellite_empty.set_startEndIndexThermalNoise(5, 10);

    // Choose median or mean to compute the thermal noise
    satellite_empty.set_thermalNoiseIsMedian(true);

    // Loop over measurements and fit each time
    double alt, lat;
    vector<double> waveform, fit_parameters_initial, fit_parameters_estimated, fit_values;
    vector<vector<double>> fit_values_matrix_obj (N, vector<double>(D, -9999));
    vector<vector<double>> fit_values_matrix_obj_hayne_skewness_2 (N, vector<double>(D, -9999));
    vector<vector<double>> fit_values_matrix_fun (N, vector<double>(D, -9999));
    vector<vector<double>> fit_param_matrix_obj (N, vector<double>(4, -9999));
    vector<vector<double>> fit_param_matrix_obj_hayne (N, vector<double>(4, -9999));

    for (unsigned int meas_idx=0; meas_idx<N; meas_idx++)
    {

        // extract the current waveform, its latitude and altitude
        waveform = waveforms[meas_idx];
        alt = alts[meas_idx][0];
        lat = lats[meas_idx][0];

        // fit_brown
        // =========
        satellite_empty.fit_brown(waveform, alt, deg2rad(lat));

        // retrievals of fit information
        fit_parameters_initial = satellite_empty.get_fit_parameters_initial();
        fit_parameters_estimated = satellite_empty.get_fit_parameters_estimated();
        fit_values = satellite_empty.get_fit_values();

        // storage for later dump in a file
        fit_values_matrix_obj[meas_idx] = fit_values;
        fit_param_matrix_obj[meas_idx] = fit_parameters_estimated;

        // fit_hayne
        // =========
        satellite_empty.fit_hayne(waveform, alt, deg2rad(lat), 2.0);

        // retrievals of fit information
        fit_parameters_initial = satellite_empty.get_fit_parameters_initial();
        fit_parameters_estimated = satellite_empty.get_fit_parameters_estimated();
        fit_values = satellite_empty.get_fit_values();

        // storage for later dump in a file
        fit_values_matrix_obj_hayne_skewness_2[meas_idx] = fit_values;
        fit_param_matrix_obj_hayne[meas_idx] = fit_parameters_estimated;

    }

    // save the data and the obj and fun fits
    string path_fit;
    string path_param;

    path_fit = "../Brown_retracker_master/data/outputs/brown_fit_object.txt";
    save_txt(fit_values_matrix_obj, path_fit);
    path_param = "../Brown_retracker_master/data/outputs/brown_param_object.txt";
    save_txt(fit_param_matrix_obj, path_param);

    path_fit = "../Brown_retracker_master/data/outputs/brown_fit_object_hayne_skewness_2.txt";
    save_txt(fit_values_matrix_obj_hayne_skewness_2, path_fit);
    path_param = "../Brown_retracker_master/data/outputs/brown_param_object_hayne_skewness_2.txt";
    save_txt(fit_param_matrix_obj, path_param);


    return 0;  // program executed nicely
}
