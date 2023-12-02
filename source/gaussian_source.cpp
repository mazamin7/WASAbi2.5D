#define _USE_MATH_DEFINES
#include <cmath>
#include "gaussian_source.h"
#include "simulation.h"


GaussianSource::GaussianSource(int x_0, int y_0, int z_0, int t_0, int points) :SoundSource(x_0, y_0, z_0, t_0, points)
{
}


GaussianSource::~GaussianSource()
{
}

double GaussianSource::SampleSpaceValue(double x, double y, double z)
{
	// 3 sigma
	double sigma = Simulation::dh_ * points() / 3;

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

	double freq_nyq = 1 / Simulation::dt_ * 0.5;
	double freq = freq_nyq * 0.05;
	double T = 1 / freq;
	double T_samples = T / Simulation::dt_;

	if ((t >= t_0()) && (t < t_0() + T_samples)) {
		double omega = 2 * M_PI * freq;

		double arg = omega * (t - t_0()) * Simulation::dt_;

		double amp = 40 * 5E0;

		val = 1 * Simulation::c0_ * Simulation::c0_ * amp * sin(arg);
	}

	return val;
}