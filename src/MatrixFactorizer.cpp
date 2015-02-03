
#include "BassoConfig.h"

#include <Eigen/Dense>
#include <fstream>
#include <string>

#include "CommandLineOptions/MatrixFactorizerOptions.hpp"
#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/VectorSetter.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "SolutionFactory/SolutionFactory.hpp"

bool solveProblem(
		Database_ptr_t &_database,
		const MatrixFactorizerOptions &_opts,
		const Eigen::MatrixXd &_matrix,
		const Eigen::MatrixXd &_rhs,
		const Eigen::VectorXd &_startingvalue,
		GeneralMinimizer::ReturnValues &_result
		)
{
	// prepare inverse problem
	InverseProblem_ptr_t inverseproblem =
			SolutionFactory::createInverseProblem(
					_opts, _matrix, _rhs);
	MinimizerFactory::instance_ptr_t minimizer =
			SolutionFactory::createMinimizer(
					_opts, inverseproblem, _database, _opts.inner_iterations);
	if (minimizer == NULL) {
		BOOST_LOG_TRIVIAL(error)
				<< "Minimizer could not be constructed, exiting.";
		return false;
	}

	// empty true solution
	SpaceElement_ptr_t truesolution =
			inverseproblem->x->getSpace()->createElement();

	// prepare start value and dual solution
	SpaceElement_ptr_t x0 = ElementCreator::create(
			inverseproblem->x->getSpace(),
			_startingvalue);
	*inverseproblem->x = x0;
	if (x0->getSpace()->getDimension() < 10)
		BOOST_LOG_TRIVIAL(debug)
			<< "Starting at x0 = " << x0;
	SpaceElement_ptr_t dualx0;
	if (_opts.dualitytype == CommandLineOptions::defaulttype) {
		dualx0 = (*inverseproblem->x->getSpace()->getDualityMapping())(x0);
	} else {
		dualx0 = inverseproblem->x->getSpace()->getDualSpace()->createElement();
		dualx0->setZero();
	}

	// and minimize
	try{
		_result = (*minimizer)(
						inverseproblem,
						x0,
						dualx0,
						truesolution);
		minimizer->resetState();
	} catch (MinimizationIllegalValue_exception &e) {
		std::cerr << "Illegal value for "
				<< *boost::get_error_info<MinimizationIllegalValue_name>(e)
				<< std::endl;
		return false;
	}

	return true;
}

void removeSmallTables(Database_ptr_t &_database)
{
	// remove tables
	_database->removeTable("per_iteration");
	_database->removeTable("overall");
	_database->removeTable("angles");
}

int main(int argc, char **argv)
{
	/// starting timing
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

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
			<< "We solve for Y=K*X with Y =\n"
			<< data << "." << std::endl;
	} else {
		BOOST_LOG_TRIVIAL(info)
					<< "We solve for Y=K*X with Y =\n"
					<< data << "." << std::endl;
	}

	// create Database
	Database_ptr_t database =
			SolutionFactory::createDatabase(opts);
	Table& loop_table = database->addTable("loop");
	Table::Tuple_t loop_tuple;
	loop_tuple.insert( std::make_pair("k", (int)data.innerSize()), Table::Parameter);
	loop_tuple.insert( std::make_pair("n", (int)data.outerSize()), Table::Parameter);
	loop_tuple.insert( std::make_pair("p", opts.normx), Table::Parameter);
	loop_tuple.insert( std::make_pair("r", opts.normy), Table::Parameter);
	loop_tuple.insert( std::make_pair("sparse_dim", (int)opts.sparse_dim), Table::Parameter);
	for (std::vector<std::string>::const_iterator iter = opts.tuple_parameters.begin();
			iter != opts.tuple_parameters.end(); iter+=2) {
		BOOST_LOG_TRIVIAL(debug)
				<< " Adding additional parameter ("
				<< *iter << "," << *(iter+1) << ") to loop tuple.";
		loop_tuple.insert( std::make_pair(*iter, *(iter+1)), Table::Parameter);
	}
	loop_tuple.insert( std::make_pair("loop_nr", (int)0), Table::Data);
	loop_tuple.insert( std::make_pair("residual", 0.), Table::Data);

	/// construct solution starting points
	Eigen::MatrixXd spectral_matrix(data.rows(), opts.sparse_dim);
	spectral_matrix.setRandom();
	if ((spectral_matrix.innerSize() > 10) || (spectral_matrix.outerSize() > 10)) {
		BOOST_LOG_TRIVIAL(trace)
				<< "Initial spectral matrix is\n" << spectral_matrix;
	} else {
		BOOST_LOG_TRIVIAL(info)
				<< "Initial spectral matrix is\n" << spectral_matrix;
	}
	Eigen::MatrixXd pixel_matrix(opts.sparse_dim, data.cols());
	pixel_matrix.setZero();

	/// iterate over the two factors
	unsigned int loop_nr = 0;
	double residual = 0.;
	bool stop_condition = loop_nr > opts.max_loops;
	while (!stop_condition) {
		// update loop count
		++loop_nr;
		loop_tuple.replace("loop_nr", (int)loop_nr);

		/// loop over pixel dimensions
		for (unsigned int pixel_dim = 0; pixel_dim < data.cols();
				++pixel_dim) {
			/// construct and solve (approximately) inverse problem
			GeneralMinimizer::ReturnValues result;
			if (!solveProblem(
					database,
					opts,
					spectral_matrix,
					data.col(pixel_dim),
					pixel_matrix.col(pixel_dim),
					result))
				return 255;
			removeSmallTables(database);
			BOOST_LOG_TRIVIAL(debug)
				<< "Resulting vector is " << *(result.m_solution);


			Eigen::VectorXd result_vector(opts.sparse_dim, 1);
			VectorSetter<Eigen::VectorXd>::set(
					result.m_solution,
					result_vector);
//			BOOST_LOG_TRIVIAL(info)
//				<< "Pixel_matrix column has dimensions "
//				<< pixel_matrix.col(pixel_dim).innerSize()
//				<< "," << pixel_matrix.col(pixel_dim).outerSize();
			pixel_matrix.col(pixel_dim) = result_vector;
		}
		// check criterion
		{
			const Eigen::MatrixXd difference_matrix =
					data - spectral_matrix * pixel_matrix;
			residual = difference_matrix.norm();
			BOOST_LOG_TRIVIAL(info)
				<< "#" << loop_nr << " 1/2, residual is " << residual;
		}

		/// loop over channel dimensions
		for (unsigned int spectral_dim = 0; spectral_dim < data.rows();
				++spectral_dim) {
			/// construct and solve (approximately) inverse problem
			GeneralMinimizer::ReturnValues result;
			if (!solveProblem(
					database,
					opts,
					pixel_matrix.transpose(),
					data.row(spectral_dim).transpose(),
					spectral_matrix.row(spectral_dim).transpose(),
					result))
				return 255;
			removeSmallTables(database);
			BOOST_LOG_TRIVIAL(debug)
				<< "Resulting vector is " << *(result.m_solution);

			Eigen::VectorXd result_vector(opts.sparse_dim,1);
			VectorSetter<Eigen::VectorXd>::set(
					result.m_solution,
					result_vector);
//			BOOST_LOG_TRIVIAL(info)
//				<< "Spectral_matrix row has dimensions "
//				<< spectral_matrix.row(spectral_dim).innerSize()
//				<< "," << spectral_matrix.row(spectral_dim).outerSize();
			spectral_matrix.row(spectral_dim) = result_vector.transpose();
		}

		// check criterion
		{
			const Eigen::MatrixXd difference_matrix =
					data - spectral_matrix * pixel_matrix;
			residual = difference_matrix.norm();
			BOOST_LOG_TRIVIAL(info)
				<< "#" << loop_nr << " 2/2, residual is " << residual;
			loop_tuple.replace("residual", residual);
			stop_condition = (residual < opts.delta)
					|| (loop_nr > opts.max_loops);
		}

		// submit loop tuple
		loop_table.addTuple(loop_tuple);
	}
	if (loop_nr > opts.max_loops)
		BOOST_LOG_TRIVIAL(error)
			<< "Maximum number of loops " << opts.max_loops
			<< " exceeded, stopping iteration.";
	else
		BOOST_LOG_TRIVIAL(info)
			<< "Loop iteration performed " << loop_nr
			<< " times.";

	/// output solution
	{
		using namespace MatrixIO;
		if (!opts.solution_factor_one_file.string().empty()) {
			std::ofstream ost(opts.solution_factor_one_file.string().c_str());
			if (ost.good())
				try {
					ost << spectral_matrix;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully write first solution factor to file.\n";
				}
			else {
				std::cerr << "Failed to open " << opts.solution_factor_one_file.string() << std::endl;
				return 255;
			}
		} else {
			std::cout << "No first solution factor file name given." << std::endl;
		}

		if (!opts.solution_factor_two_file.string().empty()) {
			std::ofstream ost(opts.solution_factor_two_file.string().c_str());
			if (ost.good())
				try {
					ost << pixel_matrix;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully write second solution factor to file.\n";
				}
			else {
				std::cerr << "Failed to open " << opts.solution_factor_two_file.string() << std::endl;
				return 255;
			}
		} else {
			std::cout << "No second solution factor file name given." << std::endl;
		}

		if (!opts.solution_product_file.string().empty()) {
			std::ofstream ost(opts.solution_product_file.string().c_str());
			if (ost.good())
				try {
					ost << spectral_matrix * pixel_matrix;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully write solution product to file.\n";
				}
			else {
				std::cerr << "Failed to open " << opts.solution_product_file.string() << std::endl;
				return 255;
			}
		} else {
			std::cout << "No solution product file name given." << std::endl;
		}

	}
	BOOST_LOG_TRIVIAL(info)
		<< "Resulting first factor transposed is\n" << spectral_matrix.transpose();
	BOOST_LOG_TRIVIAL(info)
		<< "Resulting second factor is\n" << pixel_matrix;

	const Eigen::MatrixXd product_matrix = spectral_matrix * pixel_matrix;
	if ((data.innerSize() <= 10) && (data.outerSize() <= 10)) {
		BOOST_LOG_TRIVIAL(info)
			<< "Data matrix was\n" << data;
		BOOST_LOG_TRIVIAL(info)
			<< "Product matrix is\n" << product_matrix;
		BOOST_LOG_TRIVIAL(info)
			<< "Difference matrix is\n" << data - product_matrix;
	}
	BOOST_LOG_TRIVIAL(info)
		<< "Norm of difference is " << (data - product_matrix).norm();

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	BOOST_LOG_TRIVIAL(info) << "The operation took "
			<< boost::chrono::duration<double>(timing_end - timing_start)
			<< ".";

	/// exit
	return 0;
}
