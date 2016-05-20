/*
 * GravityOptions.cpp
 *
 *  Created on: Jul 06, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include <boost/filesystem.hpp>
#include <Gravity/Options/GravityOptions.hpp>
#include <iostream>

#include "Log/Logging.hpp"

namespace po = boost::program_options;

GravityOptions::GravityOptions() :
		depth(1.),
		discretization(100)
{}

void GravityOptions::internal_init()
{
	boost::program_options::options_description desc_basso("Gravity options");

	desc_basso.add_options()
			("depth", po::value< double >(),
					"set the depth of the mass density line")
			("discretization", po::value< unsigned int >(),
					"set the number of discretization points")
			("gravity-file", po::value< boost::filesystem::path >(),
					"set the gravity field file name, i.e. parsed right-hand-side")
			("density-file", po::value< boost::filesystem::path >(),
					"set the density file, i.e. written solution")
	        ;

	desc_all.add(desc_basso);
}

void GravityOptions::internal_parse()
{
	if (vm.count("depth")) {
		depth = vm["depth"].as<double>();
		LOG(debug, "Set depth to " << depth);
	}

	if (vm.count("discretization")) {
		discretization = vm["discretization"].as<unsigned int>();
		LOG(debug, "Set number of discretization points to " << discretization);
	}

	if (vm.count("gravity-file")) {
		gravityfield_file = vm["gravity-file"].as<boost::filesystem::path>();
		LOG(debug, "Set gravity file to " << gravityfield_file);
	}

	if (vm.count("density-file")) {
		density_file = vm["density-file"].as<boost::filesystem::path>();
		LOG(debug, "Set density file to " << density_file);
	}
}

bool GravityOptions::internal_checkSensibility() const
{
	// not yet needed, we don't adhere external gravity information
//	if ((!vm.count("gravity-file")) || (!boost::filesystem::exists(gravityfield_file))) {
//		BOOST_LOG_TRIVIAL(error)
//				<< "Gravity file not set or not present.";
//		return false;
//	}

	if (!vm.count("depth")) {
		LOG(error, "Depth is not set");
		return false;
	}

	if (!vm.count("discretization")) {
		LOG(error, "Discretization is not set");
		return false;
	}

	return true;
}

void GravityOptions::internal_setSecondaryValues()
{}

void GravityOptions::internal_store(std::ostream &_output) const
{
	_output << "# [Gravity]" << std::endl;
	writeValue<double>(_output, vm,  "depth");
	writeValue<unsigned int>(_output, vm,  "discretization");
	writeValue<boost::filesystem::path>(_output, vm,  "gravity-file");
	writeValue<boost::filesystem::path>(_output, vm,  "density-file");
}
