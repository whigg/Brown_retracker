#ifndef BROWN_H
#define BROWN_H

# include <vector>
# include <map>

class Brown
{

public:

    // =======================================
    // Constructors (default and overloadings)
    // =======================================

    // default constructor (all parameters have defaults)
    Brown(int startIndexThermalNoise = 5, int endIndexThermalNoise = 10, double scaling = 1e-12,
              double timeResolution = 0.4684,
              std::vector<double> brownParamInitGuess = std::vector<double>{1.7453292519943295e-08, 65550, 60, 6.671282e-09},  // p[0] orig: 1.7453292519943295e-08, p[0] modif: 1.9132491246579563e-08, p[3] orig: 6.671282e-09.p[3] modif: 1.6678204759907602e-09
              double widthRadarPointTargetResponse = 0.513, bool ThermalNoiseComputIsMedian = true,
              double satHeight = 700000, double satLat = 0,
              double antennaBeamWidth = 0.01509709802975095);  // optim: 0.865: 0.01509709802975095; theoretical: 1.06: 0.018500490071139894, in between: 0.9: 0.015707963267948967


    // Perhaps forget about the following initialization of the object and only consider it once the program works properly

    // constructor (all parameters have defaults except the ones that could be typically changed)
    Brown(double satelliteHeight, double earthRadius);

    Brown(std::vector<double> data_points, double satelliteHeight, double earthRadius);

//    // constructor where datapoints must be submitted
//    Brown(std::vector<double> data_points);

    // constructor where dataPoints must be submitted and all options can be changed
    Brown(std::vector<double> data_points, int startIndexThermalNoise = 5, int endIndexThermalNoise = 10, double scaling = 1e-12,
              double timeResolution = 0.4684,
              std::vector<double> brownParamInitGuess = std::vector<double>{1.7453292519943295e-08, 65550, 60, 6.671282e-09},
              double widthRadarPointTargetResponse = 0.513, double thermalNoise = 0.0,
              double satelliteHeight = 700000, double earthRadius = 6371000,
              double antennaBeamWidth = 0.015707963267948967);


    // Brown(std::vector<double> data_points, );

    // Brown(double inputs);
    // get_config();


    // ==============
    // Public methods
    // ==============

    // Public get methods
    // ==================

    void print_config() const;
    int get_startIndexThermalNoise() const;
    int get_endIndexThermalNoise() const;
    int get_speedOfLight() const;
    double get_scaling() const;
    double get_timeResolution() const;
    std::vector<double> get_dataPoints() const;
    std::vector<double> get_brownParamInitGuess() const;
    double get_widthRadarPointTargetResponse() const;
    double get_thermalNoise() const;
    double get_satelliteHeight() const;
    double get_EarthRadius() const;
    // double get_skewness() const;
    double get_antennaBeamWidth() const;
    double get_gamma() const;
    double get_a() const;
    std::vector<double> get_fit_parameters_initial() const;
    std::vector<double> get_fit_parameters_estimated() const;
    std::map<std::string, double> get_config() const;
    std::vector<double> get_fit_values() const;

    // ideas of methods



    // Public set methods
    // ==================

    void set_startEndIndexThermalNoise(int startIndexThermalNoise, int endIndexThermalNoise);
    void set_thermalNoiseIsMedian(bool is_median);
    void set_brown_param_init_guess(std::vector<double> brownParamInitGuess);
    void set_antennaBeamWidth(double antennaBeamWidth);
    // double set_skewness();


    // Public main methods
    void fit_brown(std::vector<double> data, double satHeight, double satLat);
    void fit_hayne(std::vector<double> data, double satHeight, double satLat, double skewness);
    void fit(std::vector<double> data, double satHeight, double satLat, double skewness);

private:

    // ===============
    // private members
    // ===============

    unsigned int m_start;
    unsigned int m_end;
    unsigned int m_c;
    double m_scaling;
    double m_r_t;
    std::vector<double> m_data_points;
    std::vector<double> m_data_points_scaled;
    std::vector<double> m_brown_param_init_guess;
    std::vector<double> m_brown_param_estimated;
    std::map<std::string, double> m_info_fit;
    std::map<std::string, double> m_opts;
    double m_sigma_p;
    bool m_ThermalNoiseIsMedian;
    double m_T_n;
    double m_T_n_scaled;
    double m_h;
    double m_lat;
    double m_R_e;
    double m_lambda_s;
    double m_theta_0;
    double m_gamma;
    double m_a;

    // ===============
    // Private methods
    // ===============

    // Private Set methods
    // ===================
    void set_thermalNoise();
    void set_satHeight_satLat(double satHeight, double satLat);
    void set_m_a();
    void set_m_R_e();
    void set_m_sigma_p(double widthRadarPointTargetResponse);
    void set_m_scaling(double scaling);
    void set_m_r_t(double timeResolution);

    // Private get methods
    // =====================


    // Private other methods
    // ======================
    double fit_value(double t) const;





};




#endif // BROWN_H
