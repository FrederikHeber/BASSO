
#include "BassoConfig.h"

#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <boost/assign.hpp>
#include <boost/mpl/for_each.hpp>
#include <Basso/Options/BassoOptions.hpp>
#include <Basso/printCounts.hpp>

#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Elements/SpaceElementIO.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/VectorSpaceOperationCounts.hpp"
#include "Minimizations/Spaces/OperationCountMap.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriteriaFactory.hpp"
#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraints.hpp"
#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraintsFactory.hpp"
#include "Solvers/FeasibilityProblem.hpp"
#include "Solvers/InverseProblemSolver.hpp"
#include "Solvers/SolverFactory/SolverFactory.hpp"
#include "Solvers/SplitFeasibilitySolver.hpp"

using namespace boost::assign;

int main (int argc, char *argv[])
{
	// show program information
	showVersion(std::string(argv[0]));
	showCopyright();

	BassoOptions opts;
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

	// parse matrix and vector files into instances
	Eigen::VectorXd rhs;
	Eigen::VectorXd solution;
	if (!MatrixIO::parse(opts.rhs_file.string(), "rhs", rhs))
		return 255;
	const bool solution_parsed =
			MatrixIO::parse(opts.comparison_file.string(), "solution", solution);
	InverseProblem_ptr_t inverseproblem;
	switch (opts.matrix_file.size()) {
	case 1:
		{
			Eigen::MatrixXd matrix;
			if (!MatrixIO::parse(
					opts.matrix_file[0].string(),
					"matrix", matrix))
				return 255;

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
			inverseproblem = SolverFactory::createInverseProblem(
							opts, matrix, rhs);
		}
		break;
	case 2:
		{
			// parse matrix as two factors
			Eigen::MatrixXd first_factor;
			Eigen::MatrixXd second_factor;
			if (!MatrixIO::parse(
					opts.matrix_file[0].string(),
					"first factor of matrix", first_factor))
				return 255;
			if (!MatrixIO::parse(
					opts.matrix_file[1].string(),
					"second factor of matrix", second_factor))
				return 255;

			// print parsed matrix and vector if small or high verbosity requested
			if ((second_factor.innerSize() > 10) || (first_factor.outerSize() > 10)) {
				BOOST_LOG_TRIVIAL(trace)
					<< "We solve for Ax = y with A = "
					<< first_factor*second_factor << " and y = "
					<< rhs.transpose() << std::endl;
			} else {
				BOOST_LOG_TRIVIAL(info)
					<< "We solve for Ax = y with A = "
					<< first_factor*second_factor << " and y = "
					<< rhs.transpose() << std::endl;
			}

			// prepare inverse problem
			inverseproblem = SolverFactory::createInverseProblemFromFactors(
							opts, first_factor, second_factor, rhs);
		}
		break;
	}

	// set solution to zero with now known matrix dimensions
	if (!solution_parsed) {
		const unsigned int dimension = inverseproblem->A->getSourceSpace()->getDimension();
		solution = Eigen::VectorXd(dimension);
		solution.setZero();
	}

	// prepare true solution
	SpaceElement_ptr_t truesolution =
			ElementCreator::create(
					inverseproblem->x->getSpace(),
					solution);

	// create database
	Database_ptr_t database =
			SolverFactory::createDatabase(opts);

	// prepare start value and dual solution
	SpaceElement_ptr_t x0 =
			inverseproblem->x->getSpace()->createElement();
	x0->setZero();
	if (x0->getSpace()->getDimension() < 10)
		BOOST_LOG_TRIVIAL(debug)
			<< "Starting at x0 = " << x0;

	// create auxiliary constraints
	AuxiliaryConstraintsFactory constraint_factory;
	AuxiliaryConstraints::ptr_t auxiliary_constraints =
			constraint_factory.create(opts.auxiliary_constraints);

	FeasibilityProblem::ptr_t solver;
	if (auxiliary_constraints) {
		SplitFeasibilitySolver *SFP =
				new SplitFeasibilitySolver(opts);
		FeasibilityProblem::ptr_t IP(
				new InverseProblemSolver(
						inverseproblem,
						database,
						opts,
						false /* true solution calculation */)
				);
		SFP->registerFeasibilityProblem(IP);
		SFP->registerAuxiliaryConstraints(auxiliary_constraints);
		solver.reset(SFP);
	} else {
		solver.reset(
				new InverseProblemSolver(
						inverseproblem,
						database,
						opts,
						true)
				);
	}

	// and minimize
	GeneralMinimizer::ReturnValues result = (*solver)(x0, truesolution);
	if (result.status == GeneralMinimizer::ReturnValues::error)
		return 255;

	// give result
	{
		if (inverseproblem->x->getSpace()->getDimension() > 10) {
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

#ifdef USE_TIMINGS
	using namespace VectorSpaceOperations;
	// list all incurred costs
	std::cout << "===============================================" << std::endl;
	std::cout << "This minimization incurred the following costs:" << std::endl;
	{
		std::cout << "Constant time operations: "
				<< (inverseproblem->A->getSourceSpace()->getOpCounts().getTotalConstantCounts()
					+inverseproblem->A->getSourceSpace()->getDualSpace()->getOpCounts().getTotalConstantCounts()
					+inverseproblem->A->getTargetSpace()->getOpCounts().getTotalConstantCounts()
					+inverseproblem->A->getTargetSpace()->getDualSpace()->getOpCounts().getTotalConstantCounts())
				<< std::endl;
		// put all spaces into a vector
		std::vector<NormedSpace_ptr_t> spaces;
		spaces += inverseproblem->A->getSourceSpace(), inverseproblem->A->getSourceSpace()->getDualSpace();
		spaces += inverseproblem->A->getTargetSpace(), inverseproblem->A->getTargetSpace()->getDualSpace();
		printCounts DetailPrinter(spaces);
		boost::mpl::for_each<ConstantOperationCountVector_t>(DetailPrinter);
	}
	{
		std::cout << "Linear time operations: "
				<< (inverseproblem->A->getSourceSpace()->getOpCounts().getTotalLinearCounts()
					+inverseproblem->A->getSourceSpace()->getDualSpace()->getOpCounts().getTotalLinearCounts()
					+inverseproblem->A->getTargetSpace()->getOpCounts().getTotalLinearCounts()
					+inverseproblem->A->getTargetSpace()->getDualSpace()->getOpCounts().getTotalLinearCounts())
				<< std::endl;
		// put all spaces into a vector
		std::vector<NormedSpace_ptr_t> spaces;
		spaces += inverseproblem->A->getSourceSpace(), inverseproblem->A->getSourceSpace()->getDualSpace();
		spaces += inverseproblem->A->getTargetSpace(), inverseproblem->A->getTargetSpace()->getDualSpace();
		printCounts DetailPrinter(spaces);
		boost::mpl::for_each<LinearOperationCountVector_t>(DetailPrinter);
	}
	{
		std::cout << "Quadratic time operations: "
				<< (inverseproblem->A->getCount()
					+inverseproblem->A->getAdjointMapping()->getCount())
				<< std::endl;
	}
	std::cout << "===============================================" << std::endl;
#endif

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
	BOOST_LOG_TRIVIAL(info) << "The operation took "
			<< boost::chrono::duration<double>(timing_end - timing_start)
			<< ".";

	// exit
	return 0;
}
