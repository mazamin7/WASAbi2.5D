#include "sound_source.h"
#include "gaussian_source.h"
#include "simulation.h"
#include "partition.h"
#include <iostream>


SoundSource::SoundSource(int x_0, int y_0, int z_0, int points) :x_0_(x_0), y_0_(y_0), z_0_(z_0), points_(points)
{
	static int id_generator = 0;
	id_ = id_generator++;
	std::string filename;
	std::string dir_name = std::to_string(Simulation::dh_) + "_" + std::to_string(Partition::boundary_absorption_) + "_" + std::to_string(Simulation::air_absorption_);
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
		int x_0, y_0, z_0, width;
		file >> x_0 >> y_0 >> z_0 >> width;
		if (file.eof()) break;

		sources.push_back(std::make_shared<GaussianSource>(x_0 / Simulation::dh_, y_0 / Simulation::dh_, z_0 / Simulation::dh_, width / Simulation::dh_));
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
