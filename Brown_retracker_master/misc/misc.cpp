#include "misc.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>

#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
 #endif


using namespace std;

// =====================================================
double mean(std::vector<double> x, int start, int end) {

    if ((start == end) && (start == 0)) { // median of the entire array
        start = 0;
        end = x.size();
    } else {
        int size = x.size();
        // check input args
        if (start < 0) throw std::invalid_argument( "'start' cannot be negative/n" );
        if (start > end) throw std::invalid_argument( "'start' must be smaller than 'end/n'" );
        if (end > size) throw std::invalid_argument( "'end' must be smaller than 'x' length" );
        if ((end - start) < 1) throw std::invalid_argument( "end-start (i.e. subset of the vector to compute mean) must have at least length = 1./n" );
    }


    // -----computation------//
    double x_mean;
    int x_subset_size = end - start;

    // compute the sum
    double sum=0;
    for(int i = start; i < end; i++) sum += x[i];
    x_mean = sum/x_subset_size;

    return(x_mean);
}

// ======================================================
double median (std::vector<double> x, int start, int end)
{

    if ((start == end) && (start == 0)) { // median of the entire array
        start = 0;
        end = x.size();
    } else {
        int size = x.size();
        // check input args
        if (start < 0) throw std::invalid_argument( "'start' cannot be negative/n" );
        if (start > end) throw std::invalid_argument( "'start' must be smaller than 'end/n'" );
        if (end > size) throw std::invalid_argument( "'end' must be smaller than 'x' length" );
        if ((end - start) < 1) throw std::invalid_argument( "end-start (i.e. subset of the vector to compute median) must have at least length = 1./n" );
    }


    // -----computation------//
    // get interval to compute the median on
    double x_subset_median;
    int x_subset_size = end - start;

    // initialize and set size a vector containing the interval to estimate
    std::vector<double> x_subset(x_subset_size);

    // fill the vector with the intersteting sequence
    for (int i=start; i < end; i++) x_subset[i-start] = x[i];

    // sort the elements
    // TODO (only first half to speed up the process)
    std::sort(x_subset.begin(), x_subset.end());

    // compute median
    if (x_subset_size % 2 == 0) {  // for vector of even length
        x_subset_median = ( x_subset[(x_subset_size/2)-1] + x_subset[(x_subset_size/2)] ) / 2;

    } else {  // for vector of odd length
        x_subset_median = x_subset[int(std::round(x_subset_size/2))];

    }

    return x_subset_median;

}

// ======================================================
double deg2rad (double degrees) {
        return degrees * (M_PI / 180.0);
}

// ======================================================
double rad2deg (double radians) {
    return radians * (180 / M_PI);
}

// ======================================================
double EarthRadius (double lat, bool deg) {

    if (deg) lat = deg2rad(lat);  // Latitudes in radians

    const double r1 = 6378137;  // equatorial earth radius
    const double r2 = 6356752.3142;  // polar earth radius

    double part1 = pow(cos(lat), 2.0) / pow(r1, 2.0);
    double part2 = pow(sin(lat), 2.0) / pow(r2, 2.0);
    double R_e = pow(part1 + part2, -0.5);  // radius of the earth at the specific lat

    return R_e;
}




std::map<std::string, double> get_config(double *p, double *y, int n,
                                         std::vector<double> data){

    // store y in a vector container
    std::vector<double> y_vec(n);
    for (int i=0; i<n; i++) y_vec[i] = y[i];

    // initialize the constants
    double sigma_p;  // Width of the radar point target response function (unit: )
    double T_n;  // thermal noise (unit: waveform unit)
    double start, end;  // subset of the waveform to compute T_n from
    double theta_0;  // antenna beam width (unit: deg)
    const double c = 299792458;  // speed of light (m/s)
    double r_t;  // time resolution of the (mode of the) satellite
    double h;  // satellite altitude (unit: m)
    double R_e;  // earth radius
    double lambda_s;  // Skewness of the echo of an ocean surface
    double scaling;  // scaling of waveform and p[1]
    std::map<std::string, double> config;

    //------ Check the parameters ------//
    if (data.size()>0) {  // Do not use default values


        // check that right nb of optional input are given
        if (data.size() == 10) {  // right nb
            start = data[2];
            end = data[3];

        } else {  // wrong --> abort
            throw std::invalid_argument( "When 'data' argument is given, its length must be 7.\n" );
        }

        // Time resolution (unit: s)
        if (data[7] == 0) {  // default, i.e. Cryosat-2 LRM
            r_t = 0.4684/c;
        } else if (data[7] > 0) {  // other possible value (example: Cryosat-2 SAR / SARIn = 0.2342/c)
            r_t = data[7];


        } else {  // parameter does not make sens
            throw std::invalid_argument( "'data[8]' (i.e. the satellite time resolution) "
                                         "must be > 0.\n" );
        }

        // width of the radar point target response
        if (data[0] == 0) {  // MacArthur and Hayne refined solution of Brown model (used in SAMOSA3 retracker too)
            sigma_p = 0.513 * r_t;
        } else if (data[0] == 1) {  // Original Brown model constant
            sigma_p = 0.42466090014400953 * r_t;
        } else if (data[0] == 2) {  // Envisat RA-2 constant
            sigma_p = 0.53 * r_t;
        } else {  // Wrong input --> abort
            throw std::invalid_argument( "When 'data' argument is specified, 'data[0]' must be 0, 1 or 2.\n" );
        }

        // scaling of the waveform datapoints
        if (data[9] == 0) {
            scaling = 1e-12;
        } else {
            scaling = data[9];
        }



        // thermal noise method
        if (data[1] == 0) {  // median
            if ((data[2] == data[3]) && (data[2] == 0)){  // default: median computed on bins [5;10[
                start = 5;
                end = 10;
                T_n = median(y_vec, start, end);

            } else {  // median computed on bins [start;end[
                T_n = median(y_vec, start, end);
            }

        } else if (data[1] == 1) {  // mean
            if ((data[2] == data[3]) && (data[2] == 0)){  // default: bins [5;10[
                start = 5;
                end = 10;
                T_n = mean(y_vec, start, end);

            } else {  // bins [start;end[
                T_n = mean(y_vec, start, end);
            }

        } else {  // Wrong input --> abort
            throw std::invalid_argument( "When 'data' argument is specified, 'data[1]' must be 0 or 1.\n" );
        }

        // Satellite altidue
        if (data[4] == 0) {  // default
            h = 700000;

        } else if (data[4] < 0) {  // absurd value given --> abort
            throw std::invalid_argument( "'data[4]' (i.e. satellite altitude) must be > 0.\n" );

        } else {  // valid user input used
            h = data[4];
        }

        // Earth radius
        if ((data[5] >= -180) && (data[5] <= 180)) {  // ellipsoid earth computation
            R_e = EarthRadius(data[5]);

        } else {  // input does not make sens --> abort
            throw std::invalid_argument( "'data[5]' (i.e. latitude in deg) must be in [-180;180].\n" );
        }

        lambda_s = data[6];

        // Antenna beam width
        if (data[8] == 0) {  // default, i.e. mean of along and across track antenna beam width in LRM Cryosat-2 (handbook 2015: along-track=1.06deg, across-track=1.1992deg, mean=1.1296
            theta_0 = 0.9;

        } else if (data[8] > 0) {  // non default
            theta_0 = data[8];

        } else {  // absurde parameter value --> abort
            throw std::invalid_argument( "'data[8]' (i.e. the antenna beam width) must be > 0.\n " );
        }




    } else {  // Use default values
        // time resolution of the (mode of the) satellite
        r_t = 0.4684/c;
        // MacArthur and Hayne estimation of the width of the radar point target response function
        sigma_p = 0.513 * r_t;
        // scaling of waveform data points and amplitude parameter p[1]
        scaling = 1e-12;

        // therma noise computed as the median of the subset [5,10[ of the waveform
        start = 5;
        end = 10;
        T_n = median(y_vec, start, end);
        // satellite altitude
        h = 700000;
        // Earth radius
        R_e = 6371000;
        // Skewness
        lambda_s = 0;
        // Antenna beam width (deg)
        theta_0 = 1.06;

    }

    // Compute some intermediate variables independent of the model parameters and of the time delay
    double gamma = pow(sin(deg2rad(theta_0)), 2.0) / (2*log(2));
    double a = 4*c / ( gamma * h * (1 + (h / R_e)) );

    // scaling
    T_n = T_n * scaling;


    // store the config parameters
    config.insert(std::make_pair("start", start));
    config.insert(std::make_pair("end", end));
    config.insert(std::make_pair("time_resolution", r_t));
    config.insert(std::make_pair("width_radar_point_target_response", sigma_p));
    config.insert(std::make_pair("thermal_noise", T_n));
    config.insert(std::make_pair("satellite_height", h));
    config.insert(std::make_pair("earth_radius", R_e));
    config.insert(std::make_pair("skewness", lambda_s));
    config.insert(std::make_pair("antenna_beam_width", theta_0));
    config.insert(std::make_pair("gamma", gamma));
    config.insert(std::make_pair("a", a));
    config.insert(std::make_pair("scaling", scaling));

    // modifies in place the y array so that it is scaled
    for (int i=0; i < n; i++) y[i] = y[i] * scaling;

    // convert the units of the model parameters
    p[0] = deg2rad(p[0]);   // !!!!! try to devide it by 1000 to have same order of magnitude of the parameter as the others
    // (need to devide by the same scaling in the optim brown fun and plot brown fun
    p[1] = p[1]*scaling;

    p[2] = p[2] * r_t;

    return config;

}

// continue use of scaling by looking at the next functions called
// understand the issue with skewness that is commented


// ================================================================
void save_txt(std::vector<double> vector_1d, std::string path_file,
              std::vector<std::string> columns_names){

    // check the number of columns in the list
    int N_col_labels = columns_names.size();
    int N_col_vec  = 1;
    int N_row_vec = vector_1d.size();

    bool header = true;
    if (N_col_labels == 0) header = false;

    if (header) {
        if (N_col_labels != N_col_vec) {
            throw std::invalid_argument( "The number of column names given does not match "
                                         "the number of columns of the vector to be saved (should be 1 column vector)/n" );
        }
    }

    // write
    std::ofstream file_object (path_file);
    if ( (file_object.is_open()) )
    {

        // column name writing
        if (header) file_object << columns_names.at(0) << "\n";

        // data writting
        for (int i=0; i<N_row_vec; i++) file_object << vector_1d[i] << "\n";

        // closing file
        file_object.close();

        std::cout << "Vector written to output file: " << path_file << std::endl;

    }
    else std::cout << "Unable to open file: " << path_file << std::endl;

}


// ============================================================================
void save_txt(std::vector<std::vector<double>> vector_2d, std::string path_file,
              std::vector<std::string> columns_names){

    // get matrix dimensions
    int N_row_vec = vector_2d.size();
    int N_col_vec  = vector_2d[0].size();

    // check the number of columns in the list
    int N_col_labels = columns_names.size();

    bool header = true;
    if (N_col_labels == 0) header = false;

    if (header) {
        if (N_col_labels != N_col_vec) {
            throw std::invalid_argument( "The number of column names given does not match "
                                         "the number of columns of the vector to be saved (should be 2 columns vector)/n" );
        }
    }


    // write
    std::ofstream file_object (path_file);
    if ( (file_object.is_open()) )
    {

        // column names writing
        if (header) {
            file_object << columns_names.at(0);
            for (int label = 1; label < N_col_labels; label++) file_object << ";" << columns_names.at(label);
        }


        // data writting
        if (header) file_object << "\n";
        file_object  << vector_2d[0][0];
        for (int col = 1; col < N_col_vec; col++) file_object << ";" << vector_2d[0][col];

        for (int row=1; row < N_row_vec; row++) {
            file_object << "\n";
            file_object  << vector_2d[row][0];
            for (int col = 1; col < N_col_vec; col++) file_object << ";" << vector_2d[row][col];
        }

        // closing file
        file_object.close();

        std::cout << "Vector written to output file: " << path_file << std::endl;

    }
    else std::cout << "Unable to open file: " << path_file << std::endl;

}

// check file exists
bool exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}


// load txt
// ============================================================================


template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}



std::vector<std::vector<double>> load_txt(std::string path_file, char delimiter)
{

    vector<vector<double>> data;

    // check the given path exists
    if (!exists(path_file)){
        string err_message = "File (";
        err_message.append(path_file);
        err_message.append(") does not exist \n");
        throw std::invalid_argument (err_message);
    }

    // count rows and columns
    unsigned int nb_rows = count_rows(path_file);
    unsigned int nb_columns = count_columns(path_file, delimiter);

    // load the file
    data = load_file(path_file, nb_rows, nb_columns, delimiter);

    return data;
}


// Count nb of lines in file
// =========================
unsigned int count_rows(std::string path_file)
{

    std::ifstream file(path_file);
    // new lines will be skipped unless we stop it from happening:
    file.unsetf(std::ios_base::skipws);
    // count the newlines with an algorithm specialized for counting:
    unsigned int nb_rows = count(istream_iterator<char>(file), istream_iterator<char>(), '\n') + 1;


    return nb_rows;
}


// Count nb of columns in the file (assumes equal length for each row)
// ===================================================================
unsigned int count_columns(std::string path_file, char delimiter)
{

    // read first line of the file
    ifstream file;
    file.open(path_file);
    string line;

    if (file.good()) getline(file, line);

    file.close();

    // split the lines based on the seperator to then count the number of elements
    const std::string &s = line;
    vector<string> splited_row;
    splited_row = split(s, delimiter);
    unsigned int nb_columns = splited_row.size();

    return nb_columns ;

}


// load the file
// =============
std::vector<std::vector<double>> load_file(std::string path_file, unsigned int nb_rows,
                                           unsigned int nb_columns, char delimiter)
{


    std::ifstream file(path_file);
    string line;
    vector<string> splited_row;
    vector<vector<double>> data(nb_rows, vector<double>(nb_columns, -9999));

    unsigned int row = 0;
    while (std::getline(file, line))
    {
        splited_row = split(line, delimiter);
        for(unsigned int col=0; col<nb_columns; col++) data[row][col] = stod(splited_row[col]);
        row++;
    }

    return data;
}


// Get current working directory
// =============================

string get_cwd()
{
    string current_path;
    char cCurrentPath[FILENAME_MAX];

    if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
        {
            current_path = errno;
            return current_path;
        }

    cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */
    current_path = cCurrentPath;

    return current_path;

}



