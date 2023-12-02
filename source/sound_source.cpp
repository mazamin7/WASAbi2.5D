#include "sound_source.h"
#include "gaussian_source.h"
#include "simulation.h"
#include "partition.h"
#include <iostream>
#include <assert.h>


SoundSource::SoundSource(int x_0, int y_0, int z_0, int t_0, int points, double a_freq) :x_0_(x_0), y_0_(y_0), z_0_(z_0), t_0_(t_0), points_(points), a_freq_(a_freq)
{
	static int id_generator = 0;
	id_ = id_generator++;
	std::string filename;
	std::string dir_name = std::to_string(Simulation::dh_) + "_" + std::to_string(Partition::boundary_absorption_);
	filename = "./output/" + dir_name + "/source_" + std::to_string(id_) + ".txt";
	source_.open(filename, std::ios::out);
}


SoundSource::~SoundSource()
{
}

std::vector<std::shared_ptr<SoundSource>> SoundSource::ImportSources(std::string path)
{
	std::vector<std::shared_ptr<SoundSource>> sources;

	std::ifstream file;
	file.open(path, std::ifstream::in);
	while (file.good())
	{
		double x_0, y_0, z_0, width, t_0, a_freq;
		file >> x_0 >> y_0 >> z_0 >> t_0 >> width >> a_freq;
		int x_0_samples = int(x_0 / Simulation::dh_);
		int y_0_samples = int(y_0 / Simulation::dh_);
		int z_0_samples = int(z_0 / Simulation::dh_);
		int t_0_samples = int(t_0 / Simulation::dt_);
		int width_samples = int(width / Simulation::dh_);

		assert(width_samples != 0);

		if (file.eof()) break;

		sources.push_back(std::make_shared<GaussianSource>(x_0_samples, y_0_samples, z_0_samples, t_0_samples, width_samples, a_freq));
	}
	file.close();
	for (auto source : sources)
	{
		source->RecordSource();
	}
	return sources;
}

void SoundSource::RecordSource()
{
	for (int t = 0; t < Simulation::duration_ / Simulation::dt_; t++)
	{
		source_ << this->SampleTimeValue(t) << std::endl;
	}
	source_.close();
}
