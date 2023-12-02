#include "fdtd_partition.h"
#include "simulation.h"
#include <omp.h>


int FdtdPartition::GetIndex(int x, int y, int z)
{
	if (x < 0 || x >= width_)
		return width_ * height_ * depth_;
	if (y < 0 || y >= height_)
		return width_ * height_ * depth_;
	if (z < 0 || z >= depth_)
		return width_ * height_ * depth_;
	return z * height_ * width_ + y * width_ + x;
}

FdtdPartition::FdtdPartition(int xs, int ys, int zs, int w, int h, int d)
	: Partition(xs, ys, zs, w, h, d)
{
	include_self_terms_ = true;
	should_render_ = true;
	info_.type = "FDTD";

	second_order_ = true;

	int size = width_ * height_ * depth_ + 1;

	p_old_ = (double*)malloc(size * sizeof(double));
	p_ = (double*)malloc(size * sizeof(double));
	p_new_ = (double*)malloc(size * sizeof(double));

	v_ = (double*)malloc(size * sizeof(double));

	residue_ = (double*)malloc(size * sizeof(double));
	force_ = (double*)malloc(size * sizeof(double));

	memset((void*)p_old_, 0, size * sizeof(double));
	memset((void*)p_, 0, size * sizeof(double));
	memset((void*)p_new_, 0, size * sizeof(double));
	memset((void*)v_, 0, size * sizeof(double));
	memset((void*)force_, 0, size * sizeof(double));
	memset((void*)residue_, 0, size * sizeof(double));
}


FdtdPartition::~FdtdPartition()
{
	free(p_old_);
	free(p_);
	free(p_new_);
	free(v_);
	free(force_);
	free(residue_);
}

void FdtdPartition::Update_pressure()
{
	int width = width_;
	int height = height_;
	int depth = depth_;
	auto dh = Simulation::dh_;
	auto dt = Simulation::dt_;
	auto c0 = Simulation::c0_;

	// double coefs[] = { 2.0, -27.0, 270.0, -490.0, 270.0, -27.0, 2.0 };
	// auto amp = 180.0;
	double coefs[] = { 0.0, 0.0, 1.0, -2.0, 1.0, 0.0, 0.0 };
	auto amp = 1.0;

	auto alpha_abs = air_absorption_;

#pragma omp parallel for
	for (int k = 0; k < depth; k++)
	{
		//#pragma omp parallel for
		for (int j = 0; j < height; j++)
		{
			//#pragma omp parallel for
			for (int i = 0; i < width; i++)
			{
				double d2udx2 = 0.0;
				double d2udy2 = 0.0;
				double d2udz2 = 0.0;

				for (int m = 0; m < 7; m++)
				{
					d2udx2 += coefs[m] * p_[GetIndex(i + m - 3, j, k)];
					d2udy2 += coefs[m] * p_[GetIndex(i, j + m - 3, k)];
					d2udz2 += coefs[m] * p_[GetIndex(i, j, k + m - 3)];
				}

				d2udx2 /= (amp * dh * dh);
				d2udy2 /= (amp * dh * dh);
				d2udz2 /= (amp * dh * dh);


				// p_{tt}
				double term1 = 2 * p_[GetIndex(i, j, k)] - p_old_[GetIndex(i, j, k)];

				// c_0^2 Delta p
				double term3 = c0 * c0 * (d2udx2 + d2udy2 + d2udz2);	// c^2*(d2udx2+d2udy2+d2udz2)

				
				// -2 alpha p_t
				double term4 = alpha_abs * p_old_[GetIndex(i, j, k)] / dt;


				// p_{tt} = c_0^2 Delta p - 2 alpha p_t
				p_new_[GetIndex(i, j, k)] = term1 + dt * dt * (term3 + term4 + force_[GetIndex(i, j, k)]);
				p_new_[GetIndex(i, j, k)] = p_new_[GetIndex(i, j, k)] / ( 1 + dt * alpha_abs );
			}
		}
	}

	p_old_ = p_;
	p_ = p_new_;
}

void FdtdPartition::Update_velocity()
{
	return; // velocity not considered in second order
}

double* FdtdPartition::get_pressure_field()
{
	return p_;
}

double FdtdPartition::get_pressure(int x, int y, int z)
{
	return p_[GetIndex(x, y, z)];
}

void FdtdPartition::set_pressure(int x, int y, int z, double v)
{
	p_[GetIndex(x, y, z)] = v;
}

void FdtdPartition::add_to_pressure(int x, int y, int z, double v)
{
	p_[GetIndex(x, y, z)] = p_[GetIndex(x, y, z)] + v;
}

double FdtdPartition::get_velocity(int x, int y, int z)
{
	return v_[GetIndex(x, y, z)];
}

void FdtdPartition::set_velocity(int x, int y, int z, double v)
{
	v_[GetIndex(x, y, z)] = v;
}

void FdtdPartition::add_to_velocity(int x, int y, int z, double v)
{
	v_[GetIndex(x, y, z)] = v_[GetIndex(x, y, z)] + v;
}

double FdtdPartition::get_residue(int x, int y, int z)
{
	return residue_[GetIndex(x, y, z)];
}

void FdtdPartition::set_residue(int x, int y, int z, double v)
{
	residue_[GetIndex(x, y, z)] = v;
}

void FdtdPartition::add_to_residue(int x, int y, int z, double v)
{
	residue_[GetIndex(x, y, z)] = residue_[GetIndex(x, y, z)] + v;
}

double FdtdPartition::get_force(int x, int y, int z)
{
	return force_[GetIndex(x, y, z)];
}

void FdtdPartition::set_force(int x, int y, int z, double f)
{
	force_[GetIndex(x, y, z)] = f;
}

double FdtdPartition::get_force_r(int x, int y, int z)
{
	return force_[GetIndex(x, y, z)];
}

void FdtdPartition::set_force_r(int x, int y, int z, double f)
{
	force_[GetIndex(x, y, z)] = f;
}

void FdtdPartition::reset_forces()
{
	int width = width_;
	int height = height_;
	int depth = depth_;

	memset((void*)force_, 0, width * height * depth * sizeof(double));
}

void FdtdPartition::reset_residues()
{
	int width = width_;
	int height = height_;
	int depth = depth_;

	memset((void*)residue_, 0, width * height * depth * sizeof(double));
}
