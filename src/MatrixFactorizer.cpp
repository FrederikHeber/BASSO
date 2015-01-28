
#include "BassoConfig.h"

#include <Eigen/Dense>
#include <fstream>
#include <string>

#include "CommandLineOptions/MatrixFactorizerOptions.hpp"
#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "SolutionFactory/SolutionFactory.hpp"

int main(int argc, char **argv)
{
	/// some required parameters
	MatrixFactorizerOptions opts;
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

	/// parse the matrices
	Eigen::MatrixXd data;
	{
		using namespace MatrixIO;

		{
			std::ifstream ist(opts.data_file.string().c_str());
			if (ist.good())
				try {
					ist >> data;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully parse data matrix from " << opts.data_file.string() << std::endl;
					return 255;
				}
			else {
				std::cerr << "Failed to open " << opts.data_file.string() << std::endl;
				return 255;
			}

		}
	}
	// print parsed matrix and vector if small or high verbosity requested
	if ((data.innerSize() > 10) || (data.outerSize() > 10)) {
		BOOST_LOG_TRIVIAL(trace)
			<< "We solve for Y=K*X with Y = "
			<< data << "." << std::endl;
	} else {
		BOOST_LOG_TRIVIAL(info)
					<< "We solve for Y=K*X with Y = "
					<< data << "." << std::endl;
	}

	// create Database
	Database database;
	if (!opts.iteration_file.string().empty())
		database.setDatabaseFile(opts.iteration_file.string());
	database.setReplacePresentParameterTuples(opts.database_replace);

	// create Minimizer
	MinimizerFactory factory;

	/// construct solution starting points
	Eigen::MatrixXd spectral_matrix(data.rows(), opts.sparse_dim);
	Eigen::MatrixXd pixel_matrix(opts.sparse_dim, data.cols());

	/// loop over pixel dimensions
	for (unsigned int pixel_dim = 0; pixel_dim < data.cols();
			++pixel_dim) {
		/// construct and solve (approximately) inverse problem

		// prepare inverse problem
		InverseProblem_ptr_t inverseproblem =
				SolutionFactory::createInverseProblem(
						opts, spectral_matrix, data.col(pixel_dim));
		Database_ptr_t database =
				SolutionFactory::createDatabase(opts);
		MinimizerFactory::instance_ptr_t minimizer =
				SolutionFactory::createMinimizer(
						opts, inverseproblem, database, opts.inner_iterations);
		if (minimizer == NULL) {
			BOOST_LOG_TRIVIAL(error)
					<< "Minimizer could not be constructed, exiting.";
			return 255;
		}

		// empty true solution
		SpaceElement_ptr_t truesolution =
				inverseproblem->x->getSpace()->createElement();

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
	}

	/// loop over channel dimensions
	for (unsigned int spectral_dim = 0; spectral_dim < data.rows();
			++spectral_dim) {
		// construct and solve (approximately) inverse problem
	}


	/// output solution

	/// exit
	return 0;
}
