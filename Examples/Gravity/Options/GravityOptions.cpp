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
		BOOST_LOG_TRIVIAL(debug)
			<< "Set depth to " << depth;
	}

	if (vm.count("discretization")) {
		discretization = vm["discretization"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Set number of discretization points to " << discretization;
	}

	if (vm.count("gravity-file")) {
		gravityfield_file = vm["gravity-file"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Set gravity file to " << gravityfield_file;
	}

	if (vm.count("density-file")) {
		density_file = vm["density-file"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Set density file to " << density_file;
	}
}

bool GravityOptions::internal_help_conditions() const
{
	return false;
}

bool GravityOptions::internal_checkSensibility() const
{
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
