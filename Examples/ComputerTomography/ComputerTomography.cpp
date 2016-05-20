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

#include "ComputerTomography/Options/ComputerTomographyOptions.hpp"
#include "ComputerTomography/DiscretizedRadon/Backprojection.hpp"
#include "ComputerTomography/DiscretizedRadon/DiscretizedRadon.hpp"

#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Elements/SpaceElementIO.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Mappings/LinearMappingFactory.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"
#include "Minimizations/Spaces/VectorSpaceOperationCounts.hpp"
#include "Minimizations/Spaces/OperationCountMap.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriteriaFactory.hpp"
#include "Solvers/SolverFactory/SolverFactory.hpp"

using namespace boost::assign;

int main (int argc, char *argv[])
{
	// show program information
	showVersion(std::string(argv[0]));
	showCopyright();

	ComputerTomographyOptions opts;
	opts.init();

	// parse options
	try {
		opts.parse(argc, argv);
	} catch (std::exception &e) {
		std::cerr << "An error occurred: "
				<< e.what()
				<< std::endl;
		return 255;
	}


	if (opts.showHelpConditions(argv[0]))
		return 1;

	// set verbosity level
	opts.setVerbosity();

	if (!opts.checkSensibility())
		return 255;
	opts.setSecondaryValues();

	// starting timing
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	// parse vector files into instances
	Eigen::VectorXd rhs;
	Eigen::VectorXd solution;
	const double num_pixels = opts.num_pixel_x * opts.num_pixel_y;
	const double num_measurements = opts.num_angles * opts.num_offsets;
	if (!MatrixIO::parse(opts.rhs_file.string(), "rhs", rhs)) {
		rhs = Eigen::VectorXd(num_measurements);
		rhs.setZero();
	}

	if (!MatrixIO::parse(opts.comparison_file.string(), "solution", solution)) {
		solution = Eigen::VectorXd(num_pixels);
		solution.setZero();
	}

	// check that dimensions fit
	if (solution.size() != num_pixels) {
		LOG(error, "Solution has size " << solution.size()
				<< " but product of desired pixel dimensions is " << num_pixels);
		return 255;
	}
	if (rhs.size() != num_measurements) {
		LOG(error, "Right-hand-side has size " << rhs.size()
				<< " but number of measurements is " << num_measurements);
		return 255;
	}

	SpaceElement_ptr_t truesolution;

	// create the inverse problem
	InverseProblem_ptr_t inverseproblem;
	{
		NormedSpaceFactory::args_t args_SpaceX;
		NormedSpaceFactory::args_t args_SpaceY;
		if (opts.type_spacex == "lp") {
			args_SpaceX +=
					boost::any(opts.px),
					boost::any(opts.powerx);
		} else if (opts.type_spacex == "regularized_l1") {
			args_SpaceX +=
					boost::any(opts.regularization_parameter),
					boost::any(opts.powerx);
		}
		if (opts.type_spacey == "lp") {
			args_SpaceY +=
					boost::any(opts.py),
					boost::any(opts.powery);
		}
		NormedSpace_ptr_t X = NormedSpaceFactory::create(
				num_pixels, opts.type_spacex, args_SpaceX);
		NormedSpace_ptr_t Y = NormedSpaceFactory::create(
				num_measurements, opts.type_spacey, args_SpaceY);

		Mapping_ptr_t A;
		Mapping_ptr_t A_t;
		switch (opts.radon_matrix.size())
		{
		case 0:
			// internal matrix
			{
				// the discretized radon transform and backprojection
				A.reset(new DiscretizedRadon(
						X,Y,
						opts.num_pixel_x,
						opts.num_pixel_y,
						opts.num_angles,
						opts.num_offsets));
				A_t.reset(new Backprojection(
						Y->getDualSpace(),X->getDualSpace(),
						opts.num_pixel_x,
						opts.num_pixel_y,
						opts.num_angles,
						opts.num_offsets));
				// make each the adjoint of the other
				dynamic_cast<DiscretizedRadon &>(*A).setAdjointMapping(A_t);
				dynamic_cast<Backprojection &>(*A_t).setAdjointMapping(A);
				// use transpose of Backprojection as adjoint
				dynamic_cast<DiscretizedRadon &>(*A).get().setMatrix(
						dynamic_cast<Backprojection &>(*A_t).get().getMatrix().transpose());
			}
			break;
		case 1:
			// single matrix
			{
				// parse radon matrix
				Eigen::MatrixXd radon;
				if (!MatrixIO::parse(
						opts.radon_matrix[0].string(),
						"data matrix", radon))
					return 255;
				// and create a LinearMapping through the factory
				A = LinearMappingFactory::createInstance(
						X, Y, radon, false);
				A_t = A->getAdjointMapping();
			}
			break;
		case 2:
			// two matrix factors
			{
				// parse radon matrix as two factors
				Eigen::MatrixXd first_factor;
				Eigen::MatrixXd second_factor;
				if (!MatrixIO::parse(
						opts.radon_matrix[0].string(),
						"first factor of matrix", first_factor))
					return 255;
				if (!MatrixIO::parse(
						opts.radon_matrix[1].string(),
						"second factor of matrix", second_factor))
					return 255;
				// and create a TwoFactorLinearMapping through the factory
				A = LinearMappingFactory::createTwoFactorInstance(
						X, Y, first_factor, second_factor, false);
				A_t = A->getAdjointMapping();
			}
			break;
		default:
			LOG(error, "Wrong number of matrix files given.");
			return 255;
			break;
		}

		// prepare true solution
		truesolution = ElementCreator::create(X, solution);

		// then create the SpaceElement
		SpaceElement_ptr_t y;
		if (solution.isZero())
			y = ElementCreator::create(Y, rhs);
		else {
			LOG(info, "Solution given, calculating rhs from it.");
			y = (*A)(truesolution);
		}

		if (opts.noiselevel != 0.) {
			// finally, add some noise
			Eigen::VectorXd noisevector(Y->getDimension());
			// with given seed if valid
			if (opts.seed >= 0)
				srand(opts.seed);
			noisevector.setRandom();
			SpaceElement_ptr_t noise =
					ElementCreator::create(Y, noisevector);
			*noise *= y->Norm()/noise->Norm() * opts.noiselevel;
			LOG(info, "Max and min coeffs of noise are "
					<< noise->getMaxCoefficientAndIndex().first << " and " << -1.*((-1.*noise)->getMaxCoefficientAndIndex().first));
			LOG(info, "Max and min coeffs of image are "
					<< y->getMaxCoefficientAndIndex().first << " and " << -1.*((-1.*y)->getMaxCoefficientAndIndex().first));
			*y += noise;

			using namespace MatrixIO;
		}
		if (!opts.noisy_sinogram_file.string().empty()) {
			std::ofstream ost(opts.noisy_sinogram_file.string().c_str());
			if (ost.good())
				try {
					SpaceElementWriter::output(ost, y);
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully write noisy sinogram to file.\n";
				}
			else {
				std::cerr << "Failed to open " << opts.noisy_sinogram_file.string() << std::endl;
				return 255;
			}
		} else {
			std::cout << "No noisy sinogram file given." << std::endl;
		}

		inverseproblem.reset(new InverseProblem(A,X,Y,y));
	}

	// create database
	Database_ptr_t database =
			SolverFactory::createDatabase(opts);

	// create minimizer
	MinimizerFactory::instance_ptr_t minimizer =
			SolverFactory::createMinimizer(
					opts, inverseproblem, database);
	if (minimizer == NULL) {
		LOG(error, "Minimizer could not be constructed, exiting.");
		return 255;
	}

	// prepare start value and dual solution
	SpaceElement_ptr_t x0 =
			inverseproblem->x->getSpace()->createElement();
	if (x0->getSpace()->getDimension() < 10) {
		LOG(debug, "Starting at x0 = " << x0);
	}
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
		if ((inverseproblem->x->getSpace()->getDimension() > 10)) {
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
							(*inverseproblem->A)(result.m_solution));
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
	LOG(info, "The operation took "
			<< boost::chrono::duration<double>(timing_end - timing_start) << ".");

	// exit
	return 0;
}



