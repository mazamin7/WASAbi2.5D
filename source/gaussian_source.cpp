#define _USE_MATH_DEFINES
#include <cmath>
#include "gaussian_source.h"
#include "simulation.h"


GaussianSource::GaussianSource(int x_0, int y_0, int z_0, int points) :SoundSource(x_0, y_0, z_0, points)
{
}


GaussianSource::~GaussianSource()
{
}

double GaussianSource::SampleSpaceValue(double x, double y, double z)
{
	auto sigma = Simulation::dh_ * int(points() / 3);

	double arg_x = pow(M_PI * (Simulation::dh_ * (x - x_0())), 2) / (2 * sigma * sigma);
	double arg_y = pow(M_PI * (Simulation::dh_ * (y - y_0())), 2) / (2 * sigma * sigma);
	double arg_z = pow(M_PI * (Simulation::dh_ * (z - z_0())), 2) / (2 * sigma * sigma);

	auto amp = 1 / (pow(2 * M_PI, 3 / 2) * sigma * sigma * sigma);
	// amp = 1;
	auto val = amp * exp(-arg_x) * exp(-arg_y) * exp(-arg_z);
	// val = 1;
	return val;
}

double GaussianSource::SampleTimeValue(double t)
{
	double arg = pow(M_PI * ((2 * (Simulation::c0_*Simulation::dt_ / Simulation::dh_) * t) / 6 - 2.0), 2);
	return 1e9 * exp(-arg);
}
