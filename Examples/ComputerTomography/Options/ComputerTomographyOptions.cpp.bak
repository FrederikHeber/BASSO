/*
 * ComputerTomographyOptions.cpp
 *
 *  Created on: Jul 06, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include <boost/filesystem.hpp>
#include <ComputerTomography/Options/ComputerTomographyOptions.hpp>
#include <iostream>
#include <sstream>

#include "Log/Logging.hpp"

namespace po = boost::program_options;

ComputerTomographyOptions::ComputerTomographyOptions() :
		noiselevel(0.),
		num_pixel_x(0),
		num_pixel_y(0),
		num_angles(0),
		num_offsets(0),
		seed(-1)
{}

void ComputerTomographyOptions::internal_init()
{
	boost::program_options::options_description desc_basso("ComputerTomography options");

	desc_basso.add_options()
			("phantom", po::value< boost::filesystem::path >(),
					"set the file name of the solution/phantom to compare against (BregmanDistance)")
			("noise-level", po::value< double >(),
					"set noise level to disturb the projected true solution with  (relative)")
			("noisy-sinogram", po::value< boost::filesystem::path >(),
					"set the file name of the output noisy sinogram")
	        ("num-pixels-x", po::value< unsigned int >(),
	        		"set the desired number of pixels in x direction")
			("num-pixels-y", po::value< unsigned int >(),
					"set the desired number of pixels in y direction")
			("num-angles", po::value< unsigned int >(),
					"set the desired number of angle discretization steps")
			("num-offsets", po::value< unsigned int >(),
					"set the desired number of lateral offsets of detector")
			("radon-matrix", po::value< std::vector<boost::filesystem::path> >()->multitoken(),
					"set the file name(s) to parse the discretized Radon transformation matrix from, matrix factors if multiple files")
			("seed", po::value< int >(),
					"set the random number generator seed")
	        ("sinogram", po::value< boost::filesystem::path >(),
	        		"set the vector file of the right-hand side/sinogram")
			("solution", po::value< boost::filesystem::path >(),
					"set the file name to write solution vector to, i.e. x in y = A*x")
			("solution-image", po::value< boost::filesystem::path >(),
					"set the file name to write image of solution vector to, i.e. A*x")
	        ;

	desc_all.add(desc_basso);
}

void ComputerTomographyOptions::internal_parse()
{
	if (vm.count("phantom")) {
		comparison_file = vm["phantom"].as<boost::filesystem::path>();
		LOG(debug, "Parsing true solution/phantom vector from " << comparison_file);
	}

	if (vm.count("noise-level")) {
		noiselevel = vm["noise-level"].as<double>();
		LOG(debug, "Noise level set to " << noiselevel);
	}

	if (vm.count("noisy-sinogram")) {
		noisy_sinogram_file = vm["noisy-sinogram"].as<boost::filesystem::path>();
		LOG(debug, "Filename of noisy sinogram was set to " << noisy_sinogram_file);
	}

	if (vm.count("num-pixels-x")) {
		num_pixel_x = vm["num-pixels-x"].as<unsigned int>();
		LOG(debug, "Number of x pixels was set to " << num_pixel_x);
	}

	if (vm.count("num-pixels-y")) {
		num_pixel_y = vm["num-pixels-y"].as<unsigned int>();
		LOG(debug, "Number of y pixels was set to " << num_pixel_y);
	}

	if (vm.count("num-angles")) {
		num_angles = vm["num-angles"].as<unsigned int>();
		LOG(debug, "Number of angle steps was set to " << num_angles);
	}

	if (vm.count("num-offsets")) {
		num_offsets = vm["num-offsets"].as<unsigned int>();
		LOG(debug, "Number of offsets steps was set to " << num_offsets);
	}

	if (vm.count("radon-matrix")) {
		radon_matrix = vm["radon-matrix"].as< std::vector<boost::filesystem::path> >();
		std::stringstream output;
		output << radon_matrix;
		LOG(debug, "Filename(s) of radon matrix files set to " << output.str());
	}

	if (vm.count("seed")) {
		seed = vm["seed"].as<int>();
		LOG(debug, "Random number generator seed was set to " << seed);
	}

	if (vm.count("sinogram")) {
		rhs_file = vm["sinogram"].as<boost::filesystem::path>();
		LOG(debug, "Filename of sinogram was set to " << rhs_file);
	}

	if (vm.count("solution")) {
		solution_file = vm["solution"].as<boost::filesystem::path>();
		LOG(debug, "Writing solution vector to " << solution_file);
	}

	if (vm.count("solution-image")) {
		solution_image_file = vm["solution-image"].as<boost::filesystem::path>();
		LOG(debug, "Writing image of solution vector to " << solution_image_file);
	}
}

bool ComputerTomographyOptions::internal_checkSensibility() const
{
	if (!vm.count("num-pixels-x")) {
		LOG(error, "Number of pixels in x direction not set");
		return false;
	}

	if (!vm.count("num-pixels-y")) {
		LOG(error, "Number of pixels in y direction not set");
		return false;
	}

	if (!vm.count("num-angles")) {
		LOG(error, "Number of angle discretization steps not set");
		return false;
	}

	if (!vm.count("num-offsets")) {
		LOG(error, "Number of lateral offsets not set");
		return false;
	}

	if (((!vm.count("sinogram")) || (!boost::filesystem::exists(rhs_file)))
		&& ((!vm.count("phantom")) || (!boost::filesystem::exists(comparison_file)))) {
		LOG(error, "Both right-hand side and solution file not set or not present.");
		return false;

	}

	if ((vm.count("noisy-sinogram")) && (boost::filesystem::exists(noisy_sinogram_file))) {
		LOG(error, "Noisy sinogram file is already present.");
		return false;
	}

	if (!internal_checkSensibility_radonfactors())
		return false;

	return true;
}

bool ComputerTomographyOptions::internal_checkSensibility_radonfactors() const
{
	// check whether both or none are given
	switch (radon_matrix.size())
	{
	case 0:
		if (vm.count("radon-matrix")) {
			LOG(error, "Radon matrix used but not files given.");
			return false;
		}
		break;
	case 1:
		if (!boost::filesystem::exists(radon_matrix[0])) {
			LOG(error, "File of radon matrix does not exist.");
			return false;
		}
		break;
	case 2:
		if ((!boost::filesystem::exists(radon_matrix[0]))
				|| (!boost::filesystem::exists(radon_matrix[1]))) {
			LOG(error, "At least one of the radon matrix files does not exist.");
			return false;
		}
		break;
	default:
		LOG(error, "More than two Radon matrix files given.");
		return false;
		break;
	}

	return true;
}

void ComputerTomographyOptions::internal_setSecondaryValues()
{}

void ComputerTomographyOptions::internal_store(std::ostream &_output) const
{
	_output << "# [ComputerTomography]" << std::endl;
	writeValue<boost::filesystem::path>(_output, vm,  "phantom");
	writeValue<double>(_output, vm,  "noise-level");
	writeValue<boost::filesystem::path>(_output, vm,  "noisy-sinogram");
	writeValue<unsigned int>(_output, vm,  "num-pixels-x");
	writeValue<unsigned int>(_output, vm,  "num-pixels-y");
	writeValue<unsigned int>(_output, vm,  "num-angles");
	writeValue<unsigned int>(_output, vm,  "num-offsets");
	if (radon_matrix.size() != 0) {
		for (std::vector<boost::filesystem::path>::const_iterator iter = radon_matrix.begin();
				iter != radon_matrix.end();)
			_output << "\tradon-matrix = " << *(iter++) << std::endl;
	}
	writeValue<int>(_output, vm,  "seed");
	writeValue<boost::filesystem::path>(_output, vm,  "sinogram");
	writeValue<boost::filesystem::path>(_output, vm,  "solution");
	writeValue<boost::filesystem::path>(_output, vm,  "solution-image");
}
