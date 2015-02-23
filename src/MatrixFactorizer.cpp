
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
#include "Minimizations/Functions/Minimizers/MinimizerExceptions.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "SolutionFactory/SolutionFactory.hpp"

#define TRUESOLUTION 1

template <class T>
bool solveProblem(
		Database_ptr_t &_database,
		const MatrixFactorizerOptions &_opts,
		const Eigen::MatrixXd &_matrix,
		const Eigen::MatrixXd &_rhs,
		const Eigen::VectorXd &_startingvalue,
		T &_solution,
		bool _nonnegative = false
		)
{
	GeneralMinimizer::ReturnValues result;
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

#ifdef TRUESOLUTION
	// empty or true solution from diagonalization
	Eigen::MatrixXd copymatrix = _matrix;
	Eigen::JacobiSVD<Eigen::MatrixXd> svd =
			copymatrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
	const Eigen::VectorXd truesolution_vector =
			svd.solve(_rhs);
	BOOST_LOG_TRIVIAL(info)
			<< "True solution is " << truesolution_vector.transpose()
			<< " with norm "
			<< (_matrix*truesolution_vector - _rhs).norm()/_rhs.norm();
	SpaceElement_ptr_t truesolution =
			ElementCreator::create(
					inverseproblem->x->getSpace(),
					truesolution_vector);
#else 
	SpaceElement_ptr_t truesolution = 
			inverseproblem->x->getSpace()->createElement();
#endif

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
		return false;
	}
	setResultingVector(result.m_solution, _solution, _nonnegative);

	return true;
}

void removeSmallTables(Database_ptr_t &_database)
{
	// remove tables
	_database->removeTable("per_iteration");
	_database->removeTable("overall");
	_database->removeTable("angles");
}

void renormalizeMatrixByTrace(
		Eigen::MatrixXd &_matrix)
{
	const double factor = _matrix.diagonal().maxCoeff();
	if (fabs(factor) > BASSOTOLERANCE)
		_matrix *= 1./factor;
//	if (_matrix.hasNaN())
//		throw MinimizerIllegalNumber_exception()
//		<< MinimizerIllegalNumber_variablename("matrix");
}

template <class T>
void setResultingVector(
		const SpaceElement_ptr_t &_element,
		T &_vector,
		bool _nonnegative)
{
	_vector = RepresentationAdvocate::get(_element);
	if (_nonnegative)
		for (unsigned int i=0;i<_element->getSpace()->getDimension();++i)
			_vector[i] = std::max(0.,_vector[i]);
}


inline bool checkResidualCondition(
		const double _residual,
		const double _delta)
{
	return _residual < _delta;
}

inline bool checkIterationCondition(
		const unsigned int _iterations,
		const unsigned int _max_iterations)
{
	return _iterations > _max_iterations;
}

inline double calculateResidual(
		const Eigen::MatrixXd &_data,
		const Eigen::MatrixXd &_spectral_matrix,
		const Eigen::MatrixXd &_pixel_matrix
		)
{
	const Eigen::MatrixXd difference_matrix =
			_data - _spectral_matrix * _pixel_matrix;
	return difference_matrix.norm();
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
	Table& loop_overall_table = database->addTable("loop_overall");
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
	Table::Tuple_t overall_tuple = loop_tuple;
	loop_tuple.insert( std::make_pair("loop_nr", (int)0), Table::Data);
	loop_tuple.insert( std::make_pair("residual", 0.), Table::Data);
	overall_tuple.insert( std::make_pair("loops", (int)0), Table::Data);
	overall_tuple.insert( std::make_pair("residual", 0.), Table::Data);

	/// construct solution starting points
	Eigen::MatrixXd spectral_matrix(data.rows(), opts.sparse_dim);
	spectral_matrix.setRandom();
	renormalizeMatrixByTrace(spectral_matrix);
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
	double residual = calculateResidual(data, spectral_matrix, pixel_matrix);
	BOOST_LOG_TRIVIAL(info)
		<< "#" << loop_nr << " 1/2, residual is " << residual;
	loop_tuple.replace("residual", residual);
	bool stop_condition =
			checkResidualCondition(residual, opts.delta)
			|| checkIterationCondition(loop_nr, opts.max_loops);
	// submit loop tuple
	loop_table.addTuple(loop_tuple);

	while (!stop_condition) {
		// update loop count
		++loop_nr;
		loop_tuple.replace("loop_nr", (int)loop_nr);

		renormalizeMatrixByTrace(spectral_matrix);
		/// loop over pixel dimensions
		for (unsigned int pixel_dim = 0; pixel_dim < data.cols();
				++pixel_dim) {
			/// construct and solve (approximately) inverse problem
			GeneralMinimizer::ReturnValues result;
			Eigen::VectorXd pixel_matrix_col(pixel_matrix.col(pixel_dim));
			if (!solveProblem(
					database,
					opts,
					spectral_matrix,
					data.col(pixel_dim),
					pixel_matrix.col(pixel_dim),
					pixel_matrix_col,
					loop_nr >= 3))
				return 255;
			pixel_matrix.col(pixel_dim) = pixel_matrix_col;
			removeSmallTables(database);
			BOOST_LOG_TRIVIAL(info)
				<< "Resulting vector is " << pixel_matrix.col(pixel_dim).transpose();
		}
		// check criterion
		{
			residual = calculateResidual(data, spectral_matrix, pixel_matrix);
			BOOST_LOG_TRIVIAL(info)
				<< "#" << loop_nr << " 1/2, residual is " << residual;
		}

		renormalizeMatrixByTrace(pixel_matrix);
		/// loop over channel dimensions
		for (unsigned int spectral_dim = 0; spectral_dim < data.rows();
				++spectral_dim) {
			/// construct and solve (approximately) inverse problem
			GeneralMinimizer::ReturnValues result;
			Eigen::VectorXd spectral_matrix_row(spectral_matrix.row(spectral_dim));
			if (!solveProblem(
					database,
					opts,
					pixel_matrix.transpose(),
					data.row(spectral_dim).transpose(),
					spectral_matrix.row(spectral_dim).transpose(),
					spectral_matrix_row,
					loop_nr >= 3))
				return 255;
			spectral_matrix.row(spectral_dim) = spectral_matrix_row;
			removeSmallTables(database);
			BOOST_LOG_TRIVIAL(info)
				<< "Resulting vector is " << spectral_matrix.row(spectral_dim).transpose();
		}

		// check criterion
		{
			residual = calculateResidual(data, spectral_matrix, pixel_matrix);
			BOOST_LOG_TRIVIAL(info)
				<< "#" << loop_nr << " 2/2, residual is " << residual;
			loop_tuple.replace("residual", residual);
			stop_condition =
						checkResidualCondition(residual, opts.delta)
						|| checkIterationCondition(loop_nr, opts.max_loops);
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
	overall_tuple.replace("loops", (int)loop_nr);
	overall_tuple.replace("residual", residual);
	loop_overall_table.addTuple(overall_tuple);

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
