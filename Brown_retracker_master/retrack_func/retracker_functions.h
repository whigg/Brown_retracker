#ifndef RETRACKER_FUNCTIONS_H
#define RETRACKER_FUNCTIONS_H
#include <QList>
#include <map>


/**
    Compute the full Brown model as described in the 'Coastal Altimetry' book pp. 81 - 83.

    @param *p vector containing the parameters to base the optimization on.
            p[0]: offnadir mispointing (unit: radian)
            p[1]: Amplitude of the signal (unit: waveform power unit)
            p[2]: epoch (unit: s)
            p[3]: slope of the leading edge (unit: s)
    @param *y vector containing the waveform (i.e. power echo) to optimize.
    @param m number of parameters.
    @param n number of observations (i.e. the length of the waveform vector).
    @param void* config_void_ptr: void pointer (coming from a map pointer) containing information to use specific options (defaults are 0).\n
            data[0]: sigma_p estimation (width of the radar point target response function). \n
            0 = MacArthur and Heyne ; 1 = Brown ; 2 = Envisat RA-2.

            data[1]: method for T_n estimation (thermal noise). 0 = median ; 1 = mean.

            data[2]: waveform subset start (included) for T_n estimation. 0 = default = 5.

            data[3]: waveform subset end (excluded) for T_n estimation. 0 = default = 10.

            data[4]: satellite altitude (unit: m). 0 = default = 700 000.

            data[5]: satellite latitude (unit: deg) used to compute the earth radius.
            0 = default = ellipsoidal Earth radius at equator. If data[5] in [-180,180], ellipsoidal earth assumed and radius computed for the given latitude.

            data[6]: skewness of the echo of an ocean surface.
            0 = default = brown model implementation.
            Not zero = Hayne implementation with the given skewness (currently not implemented)
    @return T_n double. The value of T_n (thermal noise) is returned.
*/
void brown(double *p, double *y, int, int n, void *config_void_ptr);

double brown(double *p, double t, std::map<std::string, double> config);


#endif // RETRACKER_FUNCTIONS_H
