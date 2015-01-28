
#include "BassoConfig.h"

// A simple program that computes the square root of a number
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <boost/assign.hpp>
//#include <boost/filesystem/path.hpp>
//#include <boost/fusion/algorithm/iteration/for_each.hpp>
//#include <boost/fusion/include/for_each.hpp>
#include <boost/mpl/for_each.hpp>

#include "CommandLineOptions/BassoOptions.hpp"
#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Elements/SpaceElementWriter.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/InverseProblems/InverseProblemFactory.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/LandweberMinimizer.hpp"
#include "Minimizations/Minimizers/SequentialSubspaceMinimizer.hpp"
#include "Minimizations/Minimizers/SequentialSubspaceMinimizerNoise.hpp"
//#include "Minimizations/Minimizers/Searchspace/LastNSearchDirections.hpp"
//#include "Minimizations/Minimizers/Searchspace/SearchspaceFactory.hpp"
#include "Minimizations/Minimizers/StepWidths/DetermineStepWidthFactory.hpp"
//#include "Minimizations/Norms/Norm.hpp"
//#include "Minimizations/Norms/NormFactory.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/VectorSpaceOperationCounts.hpp"
#include "Minimizations/Spaces/OperationCountMap.hpp"

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
	InverseProblem_ptr_t inverseproblem;
	switch (opts.dualitytype) {
	case CommandLineOptions::regularizedl1norm:
		inverseproblem = InverseProblemFactory::createRegularizedL1Instance(
				opts.regularization_parameter,
				opts.powerx,
				opts.normy,
				opts.powery,
				matrix,
				rhs);
		break;
	case CommandLineOptions::defaulttype:
	default:
		inverseproblem = InverseProblemFactory::createLpInstance(
				opts.normx,
				opts.powerx,
				opts.normy,
				opts.powery,
				matrix,
				rhs);
		break;
	}

	// prepare true solution
	SpaceElement_ptr_t truesolution =
			ElementCreator::create(
					inverseproblem->x->getSpace(),
					solution);

	// prepare start value
	SpaceElement_ptr_t x0 =
			inverseproblem->x->getSpace()->createElement();
	x0->setZero();
	if (x0->getSpace()->getDimension() < 10)
		std::cout << "Starting at x0 = " << x0 << std::endl;

	// set good maximum of inner iterations
	if (opts.maxinneriter == 0)
		opts.maxinneriter = opts.N*100;

	// call minimizer
	MinimizerFactory factory;
	Database database;
	if (!opts.iteration_file.string().empty())
		database.setDatabaseFile(opts.iteration_file.string());
	database.setReplacePresentParameterTuples(opts.database_replace);
	MinimizerFactory::instance_ptr_t minimizer;
	// set regularization parameter in case of regularizedl1norm
	minimizer =
		factory.createInstance(
				opts.type,
			inverseproblem,
			opts.delta,
			opts.maxiter,
			opts.maxinneriter,
			database,
			(const enum DetermineStepWidthFactory::stepwidth_enumeration)opts.stepwidth_type,
			opts.outputsteps);
	minimizer->MaxWalltime =
			static_cast<boost::chrono::duration<double> >(opts.maxwalltime);
	minimizer->setMinLib(opts.minlib);

	// calculate initial dual solution
	SpaceElement_ptr_t dualx0 =
			(opts.dualitytype == CommandLineOptions::defaulttype) ?
			(*inverseproblem->x->getSpace()->getDualityMapping())(x0) :
			inverseproblem->x->getSpace()->getDualSpace()->createElement();

	GeneralMinimizer::ReturnValues result;
	try {
		// create instance with some specifics
		switch(opts.type) {
		case MinimizerFactory::landweber:
			static_cast<LandweberMinimizer*>(
					minimizer.get())->setC(opts.C);
			break;
		case MinimizerFactory::sequentialsubspace:
			static_cast<SequentialSubspaceMinimizer*>(minimizer.get())->setN(opts.N);
			static_cast<SequentialSubspaceMinimizer*>(
					minimizer.get())->setupdateIndexAlgorithm(
							opts.updatetype);
			static_cast<SequentialSubspaceMinimizer*>(
					minimizer.get())->setEnforceRandomMapping(
							opts.enforceRandomMapping);
			static_cast<SequentialSubspaceMinimizer*>(minimizer.get())->setInexactLinesearch(
					opts.inexactLinesearch);
			if (opts.wolfe_constants.size() == 2) {
				static_cast<SequentialSubspaceMinimizer*>(minimizer.get())->setWolfeConstants(
						opts.wolfe_constants);
				// warning in case sanity check fails
				if (!opts.inexactLinesearch)
					BOOST_LOG_TRIVIAL(warning)
						<< "Wolfe constants set although we do perform an exact line search.";
			}
			static_cast<SequentialSubspaceMinimizer*>(
					minimizer.get())->setDoCalculateAngles(
							opts.calculateAngles);
			break;
		case MinimizerFactory::sequentialsubspace_noise:
			static_cast<SequentialSubspaceMinimizerNoise*>(
					minimizer.get())->setTau(opts.tau);
			static_cast<SequentialSubspaceMinimizerNoise*>(
					minimizer.get())->setN(opts.N);
			break;
		default:
			std::cerr << "Unknown InstanceType"
				<< MinimizerFactory::getNameForType(opts.type) << "." << std::endl;
			return 255;
			break;
		}

		// and minimize
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

	// exit
	return 0;
}
