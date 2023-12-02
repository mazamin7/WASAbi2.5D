#pragma once
#include "partition.h"

class FdtdPartition : public Partition
{
	double* p_old_{ nullptr };
	double* p_{ nullptr };
	double* p_new_{ nullptr };

	double* v_{ nullptr };

	double* residue_{ nullptr };
	double* force_;

	int GetIndex(int x, int y, int z);

public:
	FdtdPartition(int xs, int ys, int zs, int w, int h, int d, double alpha_abs_);
	~FdtdPartition();

	virtual void Update_pressure();
	virtual void Update_velocity();

	virtual double* get_pressure_field();

	virtual double get_pressure(int x, int y, int z);
	virtual void set_pressure(int x, int y, int z, double v);
	virtual void add_to_pressure(int x, int y, int z, double v);

	virtual double get_velocity(int x, int y, int z);
	virtual void set_velocity(int x, int y, int z, double v);
	virtual void add_to_velocity(int x, int y, int z, double v);

	virtual double get_residue(int x, int y, int z);
	virtual void set_residue(int x, int y, int z, double v);
	virtual void add_to_residue(int x, int y, int z, double v);

	virtual double get_force(int x, int y, int z);
	virtual void set_force(int x, int y, int z, double v);

	virtual double get_force_r(int x, int y, int z);
	virtual void set_force_r(int x, int y, int z, double v);

	virtual void reset_forces();
	virtual void reset_residues();
};

