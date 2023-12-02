#pragma once
#include <vector>
#include <string>
#include <memory>
#include <fstream>

class SoundSource
{
	int id_;
	int x_0_, y_0_, z_0_, t_0_;
	int points_;
	double a_freq_;
	std::fstream source_;

public:
	SoundSource(int x_0, int y_0, int z_0, int t_0, int points, double a_freq);
	~SoundSource();

	virtual double SampleSpaceValue(double x, double y, double z) = 0;
	virtual double SampleTimeValue(double t) = 0;

	static std::vector<std::shared_ptr<SoundSource>> ImportSources(std::string path);

	void RecordSource();

	int x_0() {
		return x_0_;
	}
	int y_0() {
		return y_0_;
	}
	int z_0() {
		return z_0_;
	}
	int t_0() {
		return t_0_;
	}

	int points() {
		return points_;
	}

	double a_freq() {
		return a_freq_;
	}

	friend class Simulation;
	friend class Partition;
};

