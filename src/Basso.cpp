
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

#include "CommandLineOptions/BassoOptions.hpp"
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

struct printCounts
{
	printCounts(const std::vector<NormedSpace_ptr_t> &_spaces) :
		spaces(_spaces)
	{}

	template <class _type>
	void operator()(_type &) {
		int sum = 0;
		for (std::vector<NormedSpace_ptr_t>::const_iterator iter = spaces.begin();
				iter != spaces.end(); ++iter)
			sum += VectorSpaceOperations::getCountTiming<_type>((*iter)->getOpCounts().instance).first;
		std::cout << "\t" << VectorSpaceOperations::TypeToString<_type>()() << ": " << sum << std::endl;
	}
private:
	const std::vector<NormedSpace_ptr_t> &spaces;
};


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
	{
		using namespace MatrixIO;

		{
			std::ifstream ist(opts.matrix_file.string().c_str());
			if (ist.good())
				try {
					ist >> matrix;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully parse matrix from " << opts.matrix_file.string() << std::endl;
					return 255;
				}
			else {
				std::cerr << "Failed to open " << opts.matrix_file.string() << std::endl;
				return 255;
			}

		}
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
				solution = Eigen::VectorXd(matrix.outerSize());
				solution.setZero();
			}
		}
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
			SolutionFactory::createInverseProblem(
					opts, matrix, rhs);

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
		std::cout << "Starting at x0 = " << x0 << std::endl;
	SpaceElement_ptr_t dualx0 =
			(opts.dualitytype == CommandLineOptions::defaulttype) ?
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
