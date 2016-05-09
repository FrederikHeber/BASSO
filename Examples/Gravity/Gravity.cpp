/*
 * \file gravity.cpp
 *
 * This implements the example "gravity" from [Hansen 2010].
 *
 *  Created on: Sep 21, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <cmath>
#include <fstream>

#include <boost/chrono.hpp>
#include <Eigen/Dense>

#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "Options/GravityOptions.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/Elements/SpaceElementIO.hpp"
#include "Solvers/SolverFactory/SolverFactory.hpp"
#include "Solvers/InverseProblemSolver.hpp"

static Eigen::VectorXd getDiscretizationVector(
		const int _N)
{
	Eigen::VectorXd discretization_vector(_N);
	const double meshwidth = 1./double(_N);
	for (Eigen::VectorXd::Index i=0;i<_N;++i)
		discretization_vector(i) = (double)i/meshwidth;
	return discretization_vector;
}

Eigen::MatrixXd GravityMatrix(
		const int _N, /* target dim */
		const int _M, /* source dim */
		const double _depth
		)
{
	Eigen::MatrixXd matrix(_N,_M);
	matrix.setZero();

	Eigen::VectorXd discretization_density =
			getDiscretizationVector(_N);
	Eigen::VectorXd discretization_field =
			getDiscretizationVector(_M);

	const double sq_depth = _depth*_depth;
	const double meshwidth = 1./double(_N);
	for (Eigen::VectorXd::Index i=0;i<_M;++i) {
		const double s = discretization_field(i);
		for (Eigen::VectorXd::Index j=0;j<_N;++j) {
			const double t = discretization_density(j);
			matrix(i,j) = meshwidth * _depth/::pow(sq_depth + (s-t)*(s-t) , 3./2.);
		}
	}


	return matrix;
}

Eigen::VectorXd GravitySolution(
		const int _N
		)
{
	Eigen::VectorXd rhs(_N);

	Eigen::VectorXd::Index inner_boundary = 0.4*_N;
	for (Eigen::VectorXd::Index i=0;i<inner_boundary;++i )
		rhs(i) = 2.;
	for (Eigen::VectorXd::Index i=inner_boundary;i<_N;++i )
		rhs(i) = 1.;

	return rhs;
}

int main (int argc, char *argv[])
{
	// show program information
	showVersion(std::string(argv[0]));
	showCopyright();

	GravityOptions opts;
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

	// create database
	Database_ptr_t database =
			SolverFactory::createDatabase(opts);

	// parse matrix and vector files into instances
	const int N = opts.discretization;
	const double depth = opts.depth;
	Eigen::MatrixXd matrix = GravityMatrix(N,N, depth);
	Eigen::VectorXd solution = GravitySolution(N);
	Eigen::VectorXd rhs = matrix * solution;

	// print parsed matrix and vector if small or high verbosity requested
	if ((matrix.innerSize() > 10) || (matrix.outerSize() > 10)) {
		BOOST_LOG_TRIVIAL(trace)
			<< "We solve for Ax = y with A = "
			<< matrix << " and y = "
			<< rhs.transpose() << std::endl;
	} else {
		BOOST_LOG_TRIVIAL(info)
			<< "We solve for Ax = y with A = "
			<< matrix << " and y = "
			<< rhs.transpose() << std::endl;
	}

	// prepare inverse problem
	InverseProblem_ptr_t inverseproblem =
			SolverFactory::createInverseProblem(
					opts, matrix, rhs);

	// solving
	InverseProblemSolver solver(
			inverseproblem,
			database,
			opts,
			false /* true solution calculation */);
	SpaceElement_ptr_t solution_start =
			inverseproblem->x->getSpace()->createElement();
	GeneralMinimizer::ReturnValues result =
			solver( solution_start );
	if (result.status == GeneralMinimizer::ReturnValues::error)
		BOOST_LOG_TRIVIAL(error) << "Something went wrong during minimization.";

	// give result
	{
		if ((matrix.innerSize() > 10) || (matrix.outerSize() > 10)) {
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

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	BOOST_LOG_TRIVIAL(info) << "The operation took "
			<< boost::chrono::duration<double>(timing_end - timing_start)
			<< ".";

	// writing solution
	{
		using namespace MatrixIO;
		if (!opts.density_file.string().empty()) {
			std::ofstream ost(opts.density_file.string().c_str());
			if (ost.good())
				try {
					SpaceElementWriter::output(ost, result.m_solution);
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully write density to file.\n";
				}
			else {
				std::cerr << "Failed to open " << opts.density_file.string() << std::endl;
				return 255;
			}
		} else {
			std::cout << "No density file given." << std::endl;
		}
	}

	// exit
	return 0;
}

