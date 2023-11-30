#define _USE_MATH_DEFINES
#include <cmath>
#include "gaussian_source.h"
#include "simulation.h"
#include <assert.h>


GaussianSource::GaussianSource(int x_0, int y_0, int z_0, int t_0, int points, int samples) :SoundSource(x_0, y_0, z_0, t_0, points, samples)
{
}


GaussianSource::~GaussianSource()
{
}

double GaussianSource::SampleSpaceValue(double x, double y, double z)
{
	// 5 sigma
	double sigma = Simulation::dh_ * int(points() / 5);

	assert(sigma != 0);

	double arg_x = pow(M_PI * (Simulation::dh_ * (x - x_0())), 2) / (2 * sigma * sigma);
	double arg_y = pow(M_PI * (Simulation::dh_ * (y - y_0())), 2) / (2 * sigma * sigma);
	double arg_z = pow(M_PI * (Simulation::dh_ * (z - z_0())), 2) / (2 * sigma * sigma);

	double amp = 1 / (pow(2 * M_PI, 3 / 2) * sigma * sigma * sigma);

	double val = amp * exp(-arg_x) * exp(-arg_y) * exp(-arg_z);

	return val;
}

double GaussianSource::SampleTimeValue(double t)
{
	double val = 0.0;

	if ( (t >= t_0()) && (t < t_0() + samples())) {
		// 3 sigma
		double sigma = Simulation::dt_ * int(samples() / 3);

		double arg = pow(M_PI * (Simulation::dt_ * (t - t_0() - samples() / 2)), 2) / (2 * sigma * sigma);

		double amp = 1 / (sqrt(2 * M_PI) * sigma);

		val = 1E7 * amp * exp(-arg);
	}

	return val;
}
