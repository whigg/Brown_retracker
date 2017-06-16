# include <vector>
# include <map>
# include <stdexcept>
# include <algorithm>
# include  <cmath>
#include <fstream>
#include <iostream>

/**
    Compute the mean of a std::vector.
    'start' and 'end' defaults are 0 resulting in computing the mean on the entire vector x.

    @param x std::vector to compute the mean from.
    @param start int included starting index of the slice of x to compute the mean from.
    @param end int excluded ending index of the slice of x to compute the mean from.
    @return mean double the mean.
*/
double mean(std::vector<double> x, int start=0, int end=0);

/**
    Compute the median of a std::vector.
    'start' and 'end' defaults are 0 resulting in computing the median on the entire vector x.

    @param x std::vector to compute the median from.
    @param start int included starting index of the slice of x to compute the median from.
    @param end in excluded ending index of the slice of x to compute the median from.
    @return median double the median.
*/
double median (std::vector<double> x, int start=0, int end=10);

/**
  Converts degrees to radians

 * @brief deg2rad double
 * @param degrees double
 * @return radians double
 */
double deg2rad (double degrees);

/**
  Convverts radians to degrees

 * @brief rad2deg double
 * @param radians double
 * @return degrees double
 */
double rad2deg (double radians);


// compute the corresponding earth radius at the given latitude
/**
  Computes the radius of the Earth at a given latitude assuming and ellipsoidal Earth
 * @brief EarthRadius
 * @param lat double the latitude.
 * @param deg bool unit of 'lat'. If true, assumes degrees, else, assumes radians. default is true.
 * @return R_e double the radius the earth at a given latitude.
 */
double EarthRadius (double lat, bool deg = true);



/**
  converts the input parameters regarding methods to use, ... into parameters usable for the brown model.
  The returned map should be used in the 'brown' function as the 'config' parameter

 * @brief get_config
 * @param p  double the parameters to fit the brown model with, i.e.:
 *                  p[0]: the off-nadir mispointing angle (unit: deg)
 *                  p[1]: the amplitude of the signal (unit: power like unit, as close in amplitude to the other parameters as possible)
 *                  p[2]: epoch (unit: bins)
 *                  p[3]: slope of the leading edge (unit: s)
 * @param y  double[] the waveform power echo (same unit as p[1]). The array is modified in place so that it is appropriately scaled to ease the optimization process
 * @param n  int the number of elements in y
 * @param data  std::vector<double> the parameters to convert into useful ones to be inserted in the brown model.
 *                                  (defaults for each parameter is set to zero. If the array is not given, defaults are used for all parameters.\n
                data[0]: sigma_p estimation (width of the radar point target response function). \n
                         0 = MacArthur and Heyne ; 1 = Brown ; 2 = Envisat RA-2.
                data[1]: method for T_n estimation (thermal noise). 0 = median ; 1 = mean.
                data[2]: waveform subset start (included) for T_n estimation. 0 = default = 5.
                data[3]: waveform subset end (excluded) for T_n estimation. 0 = default = 10.
                data[4]: satellite altitude (unit: m). 0 = default = 700 000m.
                data[5]: satellite latitude (unit: deg) used to compute the earth radius.
                         If data[5] in [-180,180], ellipsoidal earth assumed and radius computed for the given latitude.
                data[6]: skewness of the echo of an ocean surface.
                         0 = default = brown model implementation.
                         Not zero = Hayne implementation with the given skewness
                data[7]: time resolution. 0 = default = cryosat-2 LRM time resolution = 0.4684/c.
                         Other given values are interpreted as is (as long as being positive numbers). Examples could be: 0.2342/c for Cryosat-2 SAR / SARIn
                data[8]: antenna beam width. 0 = default = 1.06deg = along-track antenna beam width of the CRYosat-2 LRM radar (see Handbook 2015)
                data[9]: waveform scaling (applies to waveform datapoints and amplitude of the signal (p[1]) to keep model parameters in the same range so that optimization does not fail.
                         0 = default = 1e-12. If different from 0, will be applied as the scaling.

 * @return std::map<std::string, double> the converted parameters, i.e.:
 *                  "start": starting bin used to average the thermal noise
 *                  "end": ending bin used to average the thermal noise
 *                  "time_resolution": time_resolution of the satellite (unit: s)
 *                  "width_radar_point_target_response": the width of the radar point target response function (unit: s?)
 *                  "thermal_noise": Estimated thermal noise (unit: same as y)
 *                  "satellite_height": height of the satellite when measuring the waveform (unit: m)
 *                  "earth-radius": radius of the earth
 *                  "skewness": see data[6] unchanged
 *                  "antenna_beam_width": the antenna beam width of the radar (unit: m)
 *                  "gamma": intermediate computation only needed once per waveform,
 *                           hence computed outside of the brown optim function that gets iterated over when using the levmar library.
 *                  "a": See 'gamma'

 */
std::map<std::string, double> get_config(double *p, double *y, int n, std::vector<double> data = std::vector<double>{});

/**
  Saves a 1d vector of double elements content to a txt format file comma separeted

 * @brief save_txt void
 * @param vector_1d std::vector<double>
 * @param path_file std::string
 * @param columns_names std::vector<std::string>
 * @example std::vector<double> vector_1d(100, -9999);
            std::string path_file = ".save_test.txt";
            std::vector<std::string> columns_names{"col_1"};
            save_txt(vector_1d, path_file, columns_names);
 * @return void
 */
void save_txt(std::vector<double> vector_1d, std::string path_file,
              std::vector<std::string> columns_names = std::vector<std::string>{});

/**
  Saves a 2d vector of double elements content to a txt format file comma separeted

 * @brief save_txt void
 * @param vector_2d std::vector<std::vector<double>>
 * @param path_file std::string
 * @param columns_names std::vector<std::string>
 * @example std::vector<std::vector<double>> vector_2d(100, std::vector<double> (4, -9999));
            std::string path_file = ".save_test.txt";
            std::vector<std::string> columns_names{"col_1", "col_2", "lol", "dcd"};
            save_txt(vector_2d, path_file, columns_names);
 * @return void
 */
void save_txt(std::vector<std::vector<double>> vector_2d, std::string path_file,
              std::vector<std::string> columns_names = std::vector<std::string>{});


//
template<typename Out>
void split(const std::string &s, char delim, Out result);
std::vector<std::string> split(const std::string &s, char delim);

std::vector<std::vector<double>> load_txt(std::string path_file, char delimiter);

unsigned int count_rows(std::string path_file);

unsigned int count_columns(std::string path_file, char delimiter);

std::vector<std::vector<double>> load_file(std::string path_file, unsigned int nb_rows,
                                           unsigned int nb_columns, char delimiter);

bool exists(const std::string& name);

std::string get_cwd();

