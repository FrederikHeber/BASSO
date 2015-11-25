
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
#include "Solvers/SolverFactory/SolverFactory.hpp"

using namespace boost::assign;

int main (int argc, char *argv[])
{
	// starting timing
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	BassoOptions opts;
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

	// parse matrix and vector files into instances
	Eigen::MatrixXd matrix;
	Eigen::VectorXd rhs;
	Eigen::VectorXd solution;
	if (!MatrixIO::parse(opts.matrix_file.string(), "matrix", matrix))
		return 255;
	if (!MatrixIO::parse(opts.rhs_file.string(), "rhs", rhs))
		return 255;
	if (!MatrixIO::parse(opts.comparison_file.string(), "solution", solution)) {
		solution = Eigen::VectorXd(matrix.outerSize());
		solution.setZero();
	}

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

	// prepare true solution
	SpaceElement_ptr_t truesolution =
			ElementCreator::create(
					inverseproblem->x->getSpace(),
					solution);

	// create database
	Database_ptr_t database =
			SolverFactory::createDatabase(opts);

	// create minimizer
	MinimizerFactory::instance_ptr_t minimizer =
			SolverFactory::createMinimizer(
					opts, inverseproblem, database);
	if (minimizer == NULL) {
		BOOST_LOG_TRIVIAL(error)
				<< "Minimizer could not be constructed, exiting.";
		return 255;
	}

	// prepare start value and dual solution
	SpaceElement_ptr_t x0 =
			inverseproblem->x->getSpace()->createElement();
	x0->setZero();
	if (x0->getSpace()->getDimension() < 10)
		BOOST_LOG_TRIVIAL(debug)
			<< "Starting at x0 = " << x0;
	// only for smooth spaces we may use the duality mapping
	SpaceElement_ptr_t dualx0 =
			(inverseproblem->x->getSpace()->getNorm()->isSmooth()) ?
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
