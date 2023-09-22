#include <windows.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <random>


extern "C" {
	int __declspec(dllexport) APIENTRY UserSourceDefinition(double *data);
	int __declspec(dllexport) APIENTRY UserParamNames(char *data);
}

#define PI 3.14159265358979323846

/*
For the call to UserSourceDefinition, the data is stored as follows:

ndata = the number of elements in the array data
data[0]  = the total number of values in the passed data array
data[1]  = x (to be computed by the dll and returned)
data[2]  = y (to be computed by the dll and returned)
data[3]  = z (to be computed by the dll and returned)
data[4]  = l (to be computed by the dll and returned)
data[5]  = m (to be computed by the dll and returned)
data[6]  = n (to be computed by the dll and returned)
data[7]  = relative intensity (to be computed by the dll and returned)
data[8]  = the index of the object the source is inside of

data[20] = wavelength in µm
data[21] = millimeters per unit length (1.0 for mm, 25.4 for inches, 10.0 for cm and 1000.0 for meters)
data[22] = a random number seed
data[30] = parameter 1 from user input
data[31] = parameter 2 from user input
etc... up to data[maxdata] where maxdata = int(data[0])

The DLL must compute data 1 through 7 given the other data,
and return 0 if it works; else return -1
*/

/*

Description...

*/

// Define the Bivariate Generalized Gaussian Distribution (BGGD)
double BGGD(double xx, double yy, double omega2, double alpha) {
    return std::exp(-2 * std::pow((xx * xx + yy * yy) / omega2, alpha));
}

int __declspec(dllexport) APIENTRY UserSourceDefinition(double* data) {
    // User-defined parameters
    double omega = data[30];
    double alpha = data[31];
    double rejection_factor = data[32]; // This factor controls the size of the grid for rejection sampling
    double fiber_na = data[33];

    // Set default values if input parameter is invalid
    if (omega <= 0.0)
    {
        omega = 0.1;
    }
    if (alpha < 1)
    {
        alpha = 1.0;
    }
    if (rejection_factor <= 0.0)
    {
        rejection_factor = 2.0;
    }
    if (fiber_na <= 0.0)
    {
        fiber_na = 0.39;
    }

    // Calculate omega square for efficiency
    double omega2 = omega * omega;

    // Calculate sampling interval for X and Y
    double min_val = -rejection_factor * omega;
    double max_val = rejection_factor * omega;

    // Create a random number generator engine
    std::random_device rd;
    std::mt19937 gen(rd());

    // This is a uniform distribution from which we draw X and Y
    std::uniform_real_distribution<double> draw_xy(min_val, max_val);

    // This is a uniform distribution in [0, 1] where we draw uu
    std::uniform_real_distribution<double> draw_uu(0.0, 1.0);

    double xx, yy, uu;
    bool reject;

    // Rejection sampling
    do {
        xx = draw_xy(gen);
        yy = draw_xy(gen);
        uu = draw_uu(gen);

        // Check rejection condition
        double bggd_value = BGGD(xx, yy, omega2, alpha);
        reject = (uu > bggd_value);
    } while (reject);

    // Find the half-cone angle from the fiber NA
    double half_cone_angle = std::asin(fiber_na);

    // This is a uniform distribution where we draw theta_x/y
    std::uniform_real_distribution<double> draw_theta(-half_cone_angle, half_cone_angle);

    // Generate random angles uniformly within the half-cone angle in the X-Z, and Y-Z planes
    double theta_x = draw_theta(gen);
    double theta_y = draw_theta(gen);

    // Transform theta_x/y into LMN direction cosine
    double nn = 1.0 / std::sqrt(1.0 + std::tan(theta_x) * std::tan(theta_x) + std::tan(theta_y) * std::tan(theta_y));
    double ll = nn * std::tan(theta_x);
    double mm = nn * std::tan(theta_y);

    // Return accepted ray
    data[1] = xx;
    data[2] = yy;
    data[3] = 0.0;
    data[4] = ll;
    data[5] = mm;
    data[6] = nn;
    data[7] = 1.0;

    return 0;
}

/*

For the call to UserSourceDefinition, the data is stored as follows:
data[0] is an integer indicating the number of the parameter whose name should be returned by the function.

*/
int __declspec(dllexport) APIENTRY UserParamNames(char *data)
{
	/* this function returns the name of the parameter requested */
	unsigned int ii;
	ii = (unsigned int)(unsigned char)data[0];
	strcpy_s(data, 1, "");
	if (ii == 1) strcpy_s(data, 6, "Omega");
	if (ii == 2) strcpy_s(data, 6, "Alpha");
	if (ii == 3) strcpy_s(data, 22, "Rejection grid factor");
    if (ii == 4) strcpy_s(data, 9, "Fiber NA");
	return 0;
}
