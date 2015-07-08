/*
 * ComputerTomography.cpp
 *
 *  Created on: Jul 6, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

// A simple program that computes the square root of a number
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <boost/assign.hpp>
#include <boost/mpl/for_each.hpp>

#include "ComputerTomographyBase/Options/ComputerTomographyOptions.hpp"
#include "ComputerTomographyBase/DiscretizedRadon/DiscretizedRadonMatrix.hpp"

#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Elements/SpaceElementWriter.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/VectorSpaceOperationCounts.hpp"
#include "Minimizations/Spaces/OperationCountMap.hpp"
#include "SolutionFactory/SolutionFactory.hpp"

using namespace boost::assign;

int main (int argc, char *argv[])
{
	// starting timing
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	ComputerTomographyOptions opts;
	opts.init();

	// parse options
	opts.parse(argc, argv);

	if (opts.showHelpConditions(argv[0]))
		return 1;

	// set verbosity level
	opts.setVerbosity();

	if (!opts.checkSensibility())
		return 255;
	opts.setSecondaryValues();

	// parse vector files into instances
	Eigen::VectorXd rhs;
	Eigen::VectorXd solution;
	const double num_pixels = opts.num_pixel_x * opts.num_pixel_y;
	const double num_measurements = opts.num_angles * opts.num_offsets;
	{
		using namespace MatrixIO;

		{
			std::ifstream ist(opts.rhs_file.string().c_str());
			if (ist.good())
				try {
					ist >> rhs;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully parse rhs from " << opts.rhs_file.string() << std::endl;
					return 255;
				}
			else {
				std::cerr << "Failed to open " << opts.rhs_file.string() << std::endl;
				return 255;
			}
		}

		{
			// just try to parse solution, we ignore if file does not exist
			std::ifstream ist(opts.comparison_file.string().c_str());
			if (ist.good())
				try {
					ist >> solution;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully parse solution from " << opts.comparison_file.string() << std::endl;
					return 255;
				}
			else {
				if (opts.comparison_file.string().empty())
					BOOST_LOG_TRIVIAL(debug)
							<< "No solution file was given.";
				else {
					BOOST_LOG_TRIVIAL(error)
						<< "Could not parse solution from " << opts.comparison_file.string();
					return 255;
				}
				solution = Eigen::VectorXd(num_pixels);
				solution.setZero();
			}
		}
	}
	// check that dimensions fit
	if (solution.size() != num_pixels) {
		BOOST_LOG_TRIVIAL(error)
				<< "Solution has size " << solution.size()
				<< " but product of desired pixel dimensions is "
				<< num_pixels;
		return 255;
	}
	if (rhs.size() != num_measurements) {
		BOOST_LOG_TRIVIAL(error)
				<< "Right-hand-side has size " << rhs.size()
				<< " but number of measurements is "
				<< num_measurements;
		return 255;
	}

	// setup discretized radon transform matrix
	DiscretizedRadonMatrix RadonMatrix(
			opts.num_pixel_x,
			opts.num_pixel_y,
			opts.num_angles,
			opts.num_offsets);

	// store matrix to file if desired
	{
		using namespace MatrixIO;
		if (!opts.radon_matrix.string().empty()) {
			std::ofstream ost(opts.radon_matrix.string().c_str());
			if (ost.good())
				try {
					ost << RadonMatrix.getMatrix();
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully write radon matrix to file.\n";
				}
			else {
				std::cerr << "Failed to open " << opts.radon_matrix.string() << std::endl;
				return 255;
			}
		} else {
			std::cout << "No radon matrix file given." << std::endl;
		}
	}

	// print parsed matrix and vector if small or high verbosity requested
	if ((RadonMatrix.getMatrix().innerSize() > 10) || (RadonMatrix.getMatrix().outerSize() > 10)) {
		BOOST_LOG_TRIVIAL(trace)
			<< "We solve for Ax = y with A = "
			<< RadonMatrix.getMatrix() << " and y = "
			<< rhs.transpose() << std::endl;
	} else {
		BOOST_LOG_TRIVIAL(info)
			<< "We solve for Ax = y with A = "
			<< RadonMatrix.getMatrix() << " and y = "
			<< rhs.transpose() << std::endl;
	}

	// prepare inverse problem
	InverseProblem_ptr_t inverseproblem =
			SolutionFactory::createInverseProblem(
					opts, RadonMatrix.getMatrix(), rhs);

	// prepare true solution
	SpaceElement_ptr_t truesolution =
			ElementCreator::create(
					inverseproblem->x->getSpace(),
					solution);

	// create database
	Database_ptr_t database =
			SolutionFactory::createDatabase(opts);

	// create minimizer
	MinimizerFactory::instance_ptr_t minimizer =
			SolutionFactory::createMinimizer(
					opts, inverseproblem, database, opts.maxiter);
	if (minimizer == NULL) {
		BOOST_LOG_TRIVIAL(error)
				<< "Minimizer could not be constructed, exiting.";
		return 255;
	}
	minimizer->MaxWalltime =
			static_cast<boost::chrono::duration<double> >(opts.maxwalltime);

	// prepare start value and dual solution
	SpaceElement_ptr_t x0 =
			inverseproblem->x->getSpace()->createElement();
	x0->setZero();
	if (x0->getSpace()->getDimension() < 10)
		BOOST_LOG_TRIVIAL(debug)
			<< "Starting at x0 = " << x0;
	SpaceElement_ptr_t dualx0 =
			(opts.type_spacex == "lp") ?
			(*inverseproblem->x->getSpace()->getDualityMapping())(x0) :
			inverseproblem->x->getSpace()->getDualSpace()->createElement();

	// and minimize
	GeneralMinimizer::ReturnValues result;
	try{
		result = (*minimizer)(
						inverseproblem,
						x0,
						dualx0,
						truesolution);
		minimizer->resetState();
	} catch (MinimizationIllegalValue_exception &e) {
		std::cerr << "Illegal value for "
				<< *boost::get_error_info<MinimizationIllegalValue_name>(e)
				<< std::endl;
		return 255;
	}

	// give result
	{
		if ((RadonMatrix.getMatrix().innerSize() > 10) || (RadonMatrix.getMatrix().outerSize() > 10)) {
			std::cout << "Solution after "
					<< result.NumberOuterIterations
					<< " with relative residual of " << result.residuum/inverseproblem->y->Norm()
					<< std::endl;
		} else {
			std::cout << "Solution after " << result.NumberOuterIterations
				<< " with relative residual of " << result.residuum/inverseproblem->y->Norm()
				<< " is " << std::scientific << std::setprecision(8)
				<< inverseproblem->x
				<< std::endl;
		}
	}

	// writing solution
	{
		using namespace MatrixIO;
		if (!opts.solution_file.string().empty()) {
			std::ofstream ost(opts.solution_file.string().c_str());
			if (ost.good())
				try {
					SpaceElementWriter::output(ost, result.m_solution);
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully write solution to file.\n";
				}
			else {
				std::cerr << "Failed to open " << opts.solution_file.string() << std::endl;
				return 255;
			}
		} else {
			std::cout << "No solution file given." << std::endl;
		}

		if (!opts.solution_image_file.string().empty()) {
			std::ofstream ost(opts.solution_image_file.string().c_str());
			if (ost.good())
				try {
					SpaceElementWriter::output(
							ost,
							dynamic_cast<LinearMapping &>(*inverseproblem->A.get()) * result.m_solution);
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully write image of solution to file.\n";
				}
			else {
				std::cerr << "Failed to open " << opts.solution_image_file.string() << std::endl;
				return 255;
			}
		} else {
			std::cout << "No image of solution file given." << std::endl;
		}
	}

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	BOOST_LOG_TRIVIAL(info) << "The operation took "
			<< boost::chrono::duration<double>(timing_end - timing_start)
			<< ".";

	// exit
	return 0;
}



