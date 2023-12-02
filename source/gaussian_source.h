#pragma once
#include "sound_source.h"

class GaussianSource :public SoundSource
{
public:
	GaussianSource(int x_0, int y_0, int z_0, int t_0, int points);
	~GaussianSource();

	virtual double SampleSpaceValue(double x, double y, double z);
	virtual double SampleTimeValue(double t);
};

