/*
 * MatrixFactorizer.cpp
 *
 *  Created on: Apr 7, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include <Eigen/Dense>
#include <fstream>
#include <string>

#ifdef MPI_FOUND
#include <boost/mpi/environment.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/nonblocking.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
namespace mpi = boost::mpi;

#include <deque>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#endif

#include "Options/MatrixFactorizerOptions.hpp"
#include "Database/Database.hpp"
#include "Database/Database_mock.hpp"
#include "Log/Logging.hpp"
#include "Math/Helpers.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/VectorSetter.hpp"
#include "Minimizations/Functions/Minimizers/MinimizerExceptions.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/InverseProblems/InverseProblemFactory.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/LinearMappingFactory.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"
#include "SolutionFactory/SolutionFactory.hpp"

#define TRUESOLUTION 1

#ifdef MPI_FOUND
/** Answer from
 * http://stackoverflow.com/questions/12851126/serializing-eigens-matrix-using-boost-serialization
 */
namespace boost
{
    template<class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    inline void serialize(
        Archive & ar,
        Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & t,
        const unsigned int file_version
    )
    {
        size_t rows = t.rows();
		size_t cols = t.cols();
        ar & rows;
        ar & cols;
        if( rows * cols != t.size() )
        	t.resize( rows, cols );

        ar & boost::serialization::make_array(t.data(), t.size());
    }
}


struct GlobalData
{
	GlobalData(
			const MatrixFactorizerOptions &_opts,
			const Eigen::MatrixXd &_matrix) :
		opts(_opts),
		matrix(_matrix)
	{}
	GlobalData()
	{}

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & const_cast<MatrixFactorizerOptions &>(opts);
		ar & matrix;
	}

	const MatrixFactorizerOptions opts;
	Eigen::MatrixXd matrix;
};

struct WorkPackage
{
	WorkPackage(
			const Eigen::VectorXd &_rhs,
			const Eigen::VectorXd &_solution_startvalue,
			const int _col
			) :
		rhs(_rhs),
		solution_startvalue(_solution_startvalue),
		col(_col)
	{}
	WorkPackage(
			const int _col) :
		col(_col)
	{}
	WorkPackage() :
		col(-1)
	{}
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & rhs;
		ar & solution_startvalue;
		ar & col;
	}

	Eigen::VectorXd rhs;
	Eigen::VectorXd solution_startvalue;
	int col;
};

struct WorkResult
{
	WorkResult(
			const Eigen::VectorXd &_solution,
			const bool _solve_ok,
			const int _col
			) :
		solution(_solution),
		solve_ok(_solve_ok),
		col(_col)
	{}
	WorkResult() :
		solve_ok(true),
		col(-1)
	{}

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & solution;
		ar & solve_ok;
		ar & col;
	}

	Eigen::VectorXd solution;
	bool solve_ok;
	int col;
};
#endif

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
	// SVD only gives true solution for l2 norm
	assert( inverseproblem->y->getSpace()->getNorm()->getPvalue() == 2. );

	// empty or true solution from diagonalization
	Eigen::MatrixXd copymatrix = _matrix;
	Eigen::JacobiSVD<Eigen::MatrixXd> svd =
			copymatrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
	const Eigen::VectorXd truesolution_vector =
			svd.solve(_rhs);
	BOOST_LOG_TRIVIAL(trace)
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

bool projectOntoImage(
		Database_ptr_t &_database,
		MatrixFactorizerOptions _opts,
		const Eigen::MatrixXd &_matrix,
		const Eigen::VectorXd &_rhs,
		Eigen::VectorXd &_resultingvalue
		)
{
	// use smaller delta for the projection
	_opts.delta = 1e-8;

	// require dual values
	const double dualnormx =
			Helpers::ConjugateValue(_opts.normy);
	const double dualnormy =
			Helpers::ConjugateValue(_opts.normx);
	const double dualpowerx =
			Helpers::ConjugateValue(_opts.powerx);
	const double dualpowery =
			Helpers::ConjugateValue(_opts.powery);

	// prepare right-hand side
	NormedSpace_ptr_t Ys =
			NormedSpaceFactory::createLpInstance(
					_matrix.innerSize(), dualnormy, dualpowery);
	NormedSpace_ptr_t Xs =
			NormedSpaceFactory::createLpInstance(
					_matrix.outerSize(), dualnormx, dualpowerx);
	// and the LinearMapping
	Mapping_ptr_t As =
			LinearMappingFactory::createInstance(Ys,Xs,_matrix.transpose());
	const NormedSpace &Y = *Ys->getDualSpace();

	SpaceElement_ptr_t rhs = ElementCreator::create(Y, _rhs);
	SpaceElement_ptr_t dualrhs =
			(*Y.getDualityMapping())((-1.)*rhs);
	SpaceElement_ptr_t dualmappedrhs = (-1.)*((*As)(dualrhs));

	// prepare inverse problem: Y^\ast \rightarrow X^\ast
	InverseProblem_ptr_t inverseproblem(
			new InverseProblem(As,Ys,Xs,dualmappedrhs) );

	// prepare minimizer
	const std::string algorithm_name = _opts.algorithm_name;
	const_cast<MatrixFactorizerOptions &>(_opts).algorithm_name =
			MinimizerFactory::TypeNames[MinimizerFactory::sequentialsubspace];
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
	truesolution->setZero();

	// prepare start value and dual solution
	SpaceElement_ptr_t dualy0 =
			inverseproblem->x->getSpace()->createElement();
	dualy0->setZero();
	*inverseproblem->x = dualy0;
	if (dualy0->getSpace()->getDimension() < 10)
		BOOST_LOG_TRIVIAL(debug)
			<< "Starting at dualy0 = " << dualy0;
	SpaceElement_ptr_t y0;
	if (_opts.dualitytype == CommandLineOptions::defaulttype) {
		y0 = (*inverseproblem->x->getSpace()->getDualityMapping())(dualy0);
	} else {
		y0 = inverseproblem->x->getSpace()->getDualSpace()->createElement();
		y0->setZero();
	}

	// and minimize
	{
		GeneralMinimizer::ReturnValues result;
		try{
			result = (*minimizer)(
							inverseproblem,
							dualy0,
							y0,
							truesolution);
			minimizer->resetState();
		} catch (MinimizationIllegalValue_exception &e) {
			std::cerr << "Illegal value for "
					<< *boost::get_error_info<MinimizationIllegalValue_name>(e)
					<< std::endl;
			const_cast<MatrixFactorizerOptions &>(_opts).algorithm_name = algorithm_name;
			return false;
		}
		*result.m_solution += dualrhs;
		SpaceElement_ptr_t projected_solution =
				rhs +
				(*inverseproblem->x->getSpace()->getDualityMapping())
					(result.m_solution);
		setResultingVector(projected_solution, _resultingvalue, false);
	}
	const_cast<MatrixFactorizerOptions &>(_opts).algorithm_name = algorithm_name;

	return true;
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

void renormalizeMatrixProduct(
		Eigen::MatrixXd &_matrixone,
		Eigen::MatrixXd &_matrixtwo)
{
	const double factorone = _matrixone.diagonal().maxCoeff();
	const double factortwo = _matrixtwo.diagonal().maxCoeff();
	double factor = .5*(factorone + factortwo);
	if (fabs(factor) > BASSOTOLERANCE) {
		if (factorone > factortwo)
			factor = 1./factor;
		_matrixone *= factor;
		_matrixtwo *= 1./factor;
	}
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

inline bool checkRelativeResidualCondition(
		const double _oldresidual,
		const double _residual,
		const double _delta)
{
	if (fabs(_residual) > BASSOTOLERANCE)
		return fabs(_oldresidual-_residual)/_residual < _delta;
	else
		return fabs(_oldresidual-_residual) < _delta;
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

size_t prepareParametersTable(
		Database &_database,
		const size_t _innerSize,
		const size_t _outerSize,
		const MatrixFactorizerOptions &_opts
		)
{
	Table& parameter_table = _database.addTable("parameters");
	Table::Tuple_t& parameter_tuple = parameter_table.getTuple();
	parameter_tuple.insert( std::make_pair("k", (int)_innerSize), Table::Parameter);
	parameter_tuple.insert( std::make_pair("n", (int)_outerSize), Table::Parameter);
	parameter_tuple.insert( std::make_pair("p", _opts.normx), Table::Parameter);
	parameter_tuple.insert( std::make_pair("r", _opts.normy), Table::Parameter);
	parameter_tuple.insert( std::make_pair("sparse_dim", (int)_opts.sparse_dim), Table::Parameter);
	for (std::vector<std::string>::const_iterator iter = _opts.tuple_parameters.begin();
			iter != _opts.tuple_parameters.end(); iter+=2) {
		BOOST_LOG_TRIVIAL(debug)
				<< " Adding additional parameter ("
				<< *iter << "," << *(iter+1) << ") to loop tuple.";
		parameter_tuple.insert( std::make_pair(*iter, *(iter+1)), Table::Parameter);
	}
	// we need to add it, otherwise we cannot use checks for presence
	// as they rely on complete table info
	parameter_table.addTuple(parameter_tuple);

	size_t rowid = 0;
	if (_database.isDatabaseFileGiven()) {
		// check for presence
		if (!_database.isTuplePresentInTable(parameter_table, parameter_tuple)) {
			BOOST_LOG_TRIVIAL(debug)
					<< "Parameter tuple not present, adding to table.";
			_database.writeTable(parameter_table);
		}
		// and return
		rowid = _database.getIdOfTuplePresentInTable(
					parameter_table, parameter_tuple);
		// clear table such that present tuple is not stored again
		_database.clearTable(parameter_table.getName());
		BOOST_LOG_TRIVIAL(info)
			<< "Obtaining parameter_key " << rowid;
	} else {
		// else set rowid to arbitrary value as there is no file anyway
		rowid = 1;
	}
	return rowid;
}

bool createViews(Database &_database)
{
	// write tables beforehand
	_database.writeAllTables();

	bool status = true;
	{
		// check whether tables are present and contain elements
		const Table &data_loop_table = _database.getTableConst("data_loop");
		const Table &data_loop_overall_table =
				_database.getTableConst("data_loop_overall");
		status &= !data_loop_table.empty();
		status &= !data_loop_overall_table.empty();
		// we don't  check for parameter table's non-emptiness as is
		// possibly might be if the used parameter tuple is already in
		// the database, see setParameterKey()
	}
	if (!status)
		BOOST_LOG_TRIVIAL(error)
			<< "(Some of the) Required Tables are empty, not creating views.";
	if (status) {
		std::stringstream sql;
		sql << "CREATE VIEW IF NOT EXISTS loop AS SELECT * FROM parameters p INNER JOIN data_loop d ON p.rowid = d.parameters_fk";
		BOOST_LOG_TRIVIAL(trace)
			<< "SQL: " << sql.str();
		status &= _database.executeSQLStatement(sql.str());
	}
	if (status) {
		std::stringstream sql;
		sql << "CREATE VIEW IF NOT EXISTS loop_overall AS SELECT * FROM parameters p INNER JOIN data_loop_overall d ON p.rowid = d.parameters_fk";
		BOOST_LOG_TRIVIAL(trace)
			<< "SQL: " << sql.str();
		status &= _database.executeSQLStatement(sql.str());
	}
	return status;
}

enum SolvingDirection
{
	rowwise,
	columnwise
};

 bool solve(
		const MatrixFactorizerOptions &_opts,
		const Eigen::MatrixXd &_matrix,
		const Eigen::VectorXd &_rhs,
		const Eigen::VectorXd &_solution_start,
		Eigen::VectorXd &_solution,
		const unsigned int _dim,
		const unsigned int _loop_nr
		)
{
	Database_ptr_t mock_db(new Database_mock);
	Eigen::VectorXd projected_rhs(_rhs);
	BOOST_LOG_TRIVIAL(debug)
		<< "------------------------ n=" << _dim << " ------------------";
	// project right-hand side onto range of matrix
	BOOST_LOG_TRIVIAL(trace)
			<< "Initial y_" << _dim << " is "
			<< projected_rhs.transpose();
	{
		// project y onto image of K
		if (!projectOntoImage(
				mock_db,
				_opts,
				_matrix,
				_rhs,
				projected_rhs))
			return false;
	}
	BOOST_LOG_TRIVIAL(trace)
			<< "Projected y_" << _dim << " is "
			<< projected_rhs.transpose();
	BOOST_LOG_TRIVIAL(trace)
			<< "Difference y_" << _dim << "-y'_" << _dim << " is "
			<< (_rhs-projected_rhs).transpose();
	BOOST_LOG_TRIVIAL(debug)
			<< "........................ n=" << _dim << " ..................";

	// solve inverse problem for projected right-hand side
	_solution = _solution_start;
	if (!solveProblem(
			mock_db,
			_opts,
			_matrix,
			projected_rhs,
			_solution_start,
			_solution,
			_loop_nr >= 3))
		return false;
	BOOST_LOG_TRIVIAL(trace)
		<< "Resulting vector is " << _solution.transpose();

	return true;
}

int parseOptions(
		int _argc,
		char ** _argv,
		MatrixFactorizerOptions &_opts)
{
	_opts.init();

	// parse options
	_opts.parse(_argc, _argv);

	if (_opts.showHelpConditions(_argv[0]))
		return 1;

	// set verbosity level
	_opts.setVerbosity();

	if (!_opts.checkSensibility())
		return 255;
	_opts.setSecondaryValues();

	return 0;
}

int parseDataFile(
		const std::string &_filename,
		Eigen::MatrixXd &_data)
{
	using namespace MatrixIO;

	{
		std::ifstream ist(_filename.c_str());
		if (ist.good())
			try {
				ist >> _data;
			} catch (MatrixIOStreamEnded_exception &e) {
				std::cerr << "Failed to fully parse data matrix from " << _filename << std::endl;
				return 255;
			}
		else {
			std::cerr << "Failed to open " << _filename << std::endl;
			return 255;
		}

	}

	// print parsed matrix and vector if small or high verbosity requested
	if ((_data.innerSize() > 10) || (_data.outerSize() > 10)) {
		BOOST_LOG_TRIVIAL(trace)
			<< "We solve for Y=K*X with Y =\n"
			<< _data << "." << std::endl;
	} else {
		BOOST_LOG_TRIVIAL(info)
					<< "We solve for Y=K*X with Y =\n"
					<< _data << "." << std::endl;
	}

	return 0;
}

/** This stores the connect to a database and required tables
 *
 */
struct IterationInformation
{
	/** Constructor for class IterationInformation.
	 *
	 * Initializes database, tables, and tuples.
	 *
	 * @param _opts factorization options
	 * @param _innersize inner size of matrix to factorize
	 * @param _outersize outer size of matrix to factorize
	 */
	IterationInformation(
			const MatrixFactorizerOptions &_opts,
			const unsigned int _innersize,
			const unsigned int _outersize) :
				database(SolutionFactory::createDatabase(_opts)),
				loop_table(database->addTable("data_loop")),
				loop_overall_table(database->addTable("data_loop_overall")),
				loop_tuple(loop_table.getTuple()),
				overall_tuple(loop_overall_table.getTuple())
	{
		// install parameter set and get its key
		const size_t parameter_key = prepareParametersTable(
				*database,
				_innersize,
				_outersize,
				_opts);

		loop_tuple.insert( std::make_pair("parameters_fk", (int)parameter_key), Table::Parameter);
		loop_tuple.insert( std::make_pair("loop_nr", (int)0), Table::Data);
		loop_tuple.insert( std::make_pair("residual", 0.), Table::Data);

		overall_tuple.insert( std::make_pair("parameters_fk", (int)parameter_key), Table::Parameter);
		overall_tuple.insert( std::make_pair("loops", (int)0), Table::Data);
		overall_tuple.insert( std::make_pair("residual", 0.), Table::Data);
	}

	/** Destructor for class IterationInformation.
	 *
	 * This creates views after (some) data has been written to.
	 *
	 */
	~IterationInformation()
	{
		// add views after tables have been created and filled
		if (!createViews(*database))
			BOOST_LOG_TRIVIAL(warning)
				<< "Could not create overall or per_iteration views.";
	}

	//!> enumeration of all tables contained in this struct
	enum TableType {
		LoopTable,
		OverallTable
	};

	/** Replaces a value of key \a _keyname in the tuple of the respective
	 * matrix \a _table.
	 *
	 * \param _type type of table to add to.
	 * \param _keyname name of key whose value to replace
	 * \param _value value to replace
	 */
	template <typename T>
	void replace(
			const enum TableType _table,
			const std::string &_keyname,
			const T _value
			)
	{
		switch(_table) {
		case LoopTable:
			loop_tuple.replace(_keyname, _value);
			break;
		case OverallTable:
			overall_tuple.replace(_keyname, _value);
			break;
		default:
			BOOST_LOG_TRIVIAL(error)
				<< "Unknown table in IterationInformation::replace()";
		}
	}

	/** Adds the currently present data tuple as it is for the respective
	 * matrix \a _table to the table.
	 *
	 * \param _type type of table to add to.
	 */
	void addTuple(
			const enum TableType _table
			)
	{
		switch(_table) {
		case LoopTable:
			loop_table.addTuple(loop_tuple);
			break;
		case OverallTable:
			loop_overall_table.addTuple(overall_tuple);
			break;
		default:
			BOOST_LOG_TRIVIAL(error)
				<< "Unknown table in IterationInformation::addTuple()";
		}
	}

private:
	//!> database connection
	Database_ptr_t database;
	//!> table with per loop information
	Table& loop_table;
	//!> table with overall minimization information
	Table& loop_overall_table;
	//!> data tuple for loop_table
	Table::Tuple_t& loop_tuple;
	//!> data tuple for loop_overall_table
	Table::Tuple_t& overall_tuple;
};

void constructStartingMatrices(
		Eigen::MatrixXd &_spectral_matrix,
		Eigen::MatrixXd &_pixel_matrix
		)
{
	_spectral_matrix.setRandom();
	renormalizeMatrixByTrace(_spectral_matrix);
	if ((_spectral_matrix.innerSize() > 10) || (_spectral_matrix.outerSize() > 10)) {
		BOOST_LOG_TRIVIAL(trace)
				<< "Initial spectral matrix is\n" << _spectral_matrix;
	} else {
		BOOST_LOG_TRIVIAL(info)
				<< "Initial spectral matrix is\n" << _spectral_matrix;
	}
	_pixel_matrix.setZero();
}

#ifdef MPI_FOUND

enum MPITags {
	InitialId,
	ColumnWise,
	MAX_MPITags
};

void sendTerminate(
		mpi::communicator &world)
{
	// send full terminate signal
	// send round global information
	BOOST_LOG_TRIVIAL(debug)
			<< "#0 - broadcasting full terminate.";
	bool full_terminate = true;
	mpi::broadcast(world, full_terminate, 0);
}

bool handleResult(
		const mpi::status &_result,
		const std::vector<WorkResult> &_results,
		Eigen::MatrixXd &_solution
		)

{
	bool solve_ok = true;
	if (_result.error() != 0) {
		// look at stop_condition and store solution
		const int id = _result.source();
		solve_ok &= _results[id].solve_ok;
		_solution.col(_results[id].col) = _results[id].solution;
		if ((_solution.innerSize() > 10) || (_solution.outerSize() > 10)) {
			BOOST_LOG_TRIVIAL(trace)
					<< "Got solution for col #" << _results[id].col
					<< "\n" << _results[id].solution.transpose();
		} else {
			BOOST_LOG_TRIVIAL(debug)
					<< "Got solution for col #" << _results[id].col
					<< "\n" << _results[id].solution.transpose();
		}
	} else
		solve_ok = false;

	return solve_ok;
}

bool solveMaster(
		mpi::communicator &world,
		const MatrixFactorizerOptions &_opts,
		const Eigen::MatrixXd &_matrix,
		const Eigen::MatrixXd &_rhs,
		Eigen::MatrixXd &_solution,
		const unsigned int _loop_nr
		)
{
	// send round that we don't terminate yet
	BOOST_LOG_TRIVIAL(debug)
			<< "#0 - broadcasting no terminate.";
	bool full_terminate = false;
	mpi::broadcast(world, full_terminate, 0);

	// gather all available workers
	std::deque<int> AvailableWorkers;
	for (int i=1;i<world.size();++i) {
		int id = -1;
		BOOST_LOG_TRIVIAL(debug)
				<< "#0 - waiting for id " << i << ".";
		world.recv(i, InitialId, id);
		assert( id != -1 );
		AvailableWorkers.push_back(id);
	}

	// send round global information
	BOOST_LOG_TRIVIAL(debug)
			<< "#0 - broadcasting global data.";
	GlobalData global_data(_opts, _matrix);
	mpi::broadcast(world, global_data, 0);
	unsigned int loop_nr = _loop_nr;
	mpi::broadcast(world, loop_nr, 0);

	// allocate variables for non-blocking communication
	std::vector<WorkPackage> packages(world.size());
	std::vector<WorkResult> results(world.size());
	std::deque<mpi::request> uncompleted_requests;

	// go through all columns of rhs
	bool continue_condition = true;
	for (int col = 0;
			(continue_condition) && (col < _rhs.cols());
			++col) {
		// send current column of rhs as work package
		const size_t freeworker = AvailableWorkers.front();
		AvailableWorkers.pop_front();
		packages[freeworker] =
				WorkPackage(_rhs.col(col), _solution.col(col), col);
		BOOST_LOG_TRIVIAL(debug)
				<< "#0 - sending work package " << col
				<< " to " << freeworker << ".";
		world.isend(freeworker, ColumnWise, packages[freeworker]);

		// receive solution in a non-blocking manner
		BOOST_LOG_TRIVIAL(debug)
				<< "#0 - launching receiving work package " << col
				<< " from " << freeworker << ".";
		mpi::request request_object =
				world.irecv(freeworker, ColumnWise, results[freeworker]);
		uncompleted_requests.push_back(request_object);

		// wait for free workers if non available
		if ((AvailableWorkers.empty()) && (!uncompleted_requests.empty())) {
			// wait for any result
			BOOST_LOG_TRIVIAL(debug)
					<< "#0 - waiting for any work result.";
			std::pair<
				mpi::status,
				std::deque<mpi::request>::iterator > result =
					mpi::wait_any(
							uncompleted_requests.begin(),
							uncompleted_requests.end());
			uncompleted_requests.erase(result.second);
			BOOST_LOG_TRIVIAL(debug)
					<< "#0 - received work result from "
					<< result.first.source() << ".";
			// place worker back into available deque
			AvailableWorkers.push_back(result.first.source());
			// store solution
			continue_condition &=
					handleResult(result.first, results, _solution);
		}
	}
	// wait for all remaining results and handle incoming results
	BOOST_LOG_TRIVIAL(debug)
			<< "#0 - waiting for all remaining work results.";
	{
		std::vector<mpi::status> completed_requests;
		mpi::wait_all(
				uncompleted_requests.begin(),
				uncompleted_requests.end(),
				std::back_inserter(completed_requests));
		uncompleted_requests.clear();
		for (std::vector<mpi::status>::const_iterator iter = completed_requests.begin();
				iter != completed_requests.end(); ++iter) {
			// store solution
			continue_condition &=
					handleResult(*iter, results, _solution);
		}
	}

	// send terminate signal by empty work package
	for (int i=1;i<world.size();++i) {
		BOOST_LOG_TRIVIAL(debug)
				<< "#0 - sending termination signal to " << i << ".";
		world.send(i, ColumnWise, WorkPackage());
	}

	return continue_condition;
}

/** Worker function for arbitrary "slave" process to work on inverse
 * problems for solving matrix factorization problem.
 *
 * Workers operate in two nested loops:
 * -# the inner loop is for handling the columns of a single expression
 *    \f$ Y_i = K*X_i \f$ or $\f$ Y^t_i = X^t * K^t_i \f$, respectively.
 * -# the outer loop is for the alternating between the two matrix factors
 *    K and X. It resets the inner loop with respect to the data that is
 *    "global" therein.
 *
 * @param _world mpi communicator to receive work packages and send results
 */
void worker(mpi::communicator &world)
{
	bool full_terminate = false;
	while (!full_terminate) {
		// check whether we have to terminate
		mpi::broadcast(world, full_terminate, 0);
		if (full_terminate)
			break;

		// send id
		BOOST_LOG_TRIVIAL(debug)
				<< "#" << world.rank() << " - initiating by sending id.";
		world.send(0,InitialId,world.rank());

		// get global information
		BOOST_LOG_TRIVIAL(debug)
				<< "#" << world.rank() << " - getting global data.";
		GlobalData global_data;
		mpi::broadcast(world, global_data, 0);
		unsigned int loop_nr = 0;
		mpi::broadcast(world, loop_nr, 0);
		const Eigen::MatrixXd &matrix = global_data.matrix;
		const MatrixFactorizerOptions &opts = global_data.opts;

		if ((matrix.innerSize() > 10) || (matrix.outerSize() > 10)) {
			BOOST_LOG_TRIVIAL(trace)
					<< "#" << world.rank() << " - got global data\n"
					<< matrix;
		} else {
			BOOST_LOG_TRIVIAL(debug)
					<< "#" << world.rank() << " - got global data\n"
					<< matrix;
		}

		// we stop working only when we get the termination signal
		// from the master
		bool terminate = false;
		while ((!terminate) && (!full_terminate)) {
			/// wait for receiving data
			WorkPackage package;
			BOOST_LOG_TRIVIAL(debug)
					<< "#" << world.rank() << " - receiving next work package.";
			world.recv(0, ColumnWise, package);
			const int col = package.col;
			terminate |= (col == -1);
			if (!terminate) {
				const Eigen::VectorXd &rhs = package.rhs;
				const Eigen::VectorXd &solution_startvalue =
						package.solution_startvalue;

				if ((matrix.innerSize() > 10) || (matrix.outerSize() > 10)) {
					BOOST_LOG_TRIVIAL(trace)
							<< "#" << world.rank() << " got problem rhs, col "
							<< col << "\n" << rhs.transpose();
				} else {
					BOOST_LOG_TRIVIAL(debug)
							<< "#" << world.rank() << " got problem rhs, col "
							<< col << "\n" << rhs.transpose();
				}

				/// work on data
				BOOST_LOG_TRIVIAL(debug)
						<< "#" << world.rank() << " - working.";
				Eigen::VectorXd solution;
				const bool solve_ok =
						solve(opts,
								matrix,
								rhs,
								solution_startvalue,
								solution,
								col,
								loop_nr);

				if ((matrix.innerSize() > 10) || (matrix.outerSize() > 10)) {
					BOOST_LOG_TRIVIAL(trace)
							<< "#" << world.rank() << " sending solution, col "
							<< col << "\n" << solution.transpose();
				} else {
					BOOST_LOG_TRIVIAL(debug)
							<< "#" << world.rank() << " sending solution, col "
							<< col << "\n" << solution.transpose();
				}

				/// return result
				BOOST_LOG_TRIVIAL(debug)
						<< "#" << world.rank() << " - sending solution.";
				WorkResult result(solution, solve_ok, col);
				world.send(0, ColumnWise, result);
			} else {
				BOOST_LOG_TRIVIAL(debug)
						<< "#" << world.rank() << " - terminating.";
			}
		}
	}
}

#endif

int MatrixFactorization(
		const MatrixFactorizerOptions &_opts,
		const Eigen::MatrixXd &_data,
		IterationInformation &_info
#ifdef MPI_FOUND
		, mpi::communicator &world
#endif
		)
{
	/// construct solution starting points
	Eigen::MatrixXd spectral_matrix(_data.rows(), _opts.sparse_dim);
	Eigen::MatrixXd pixel_matrix(_opts.sparse_dim, _data.cols());
	constructStartingMatrices(spectral_matrix, pixel_matrix);

	/// iterate over the two factors
	unsigned int loop_nr = 0;
	double old_residual = 0.;
	double residual = calculateResidual(_data, spectral_matrix, pixel_matrix);
	BOOST_LOG_TRIVIAL(info)
		<< "#" << loop_nr << " 1/2, residual is " << residual;
	_info.replace(IterationInformation::LoopTable, "residual", residual);
	bool stop_condition =
			checkRelativeResidualCondition(old_residual, residual, _opts.residual_threshold)
			|| checkResidualCondition(residual, _opts.residual_threshold)
			|| checkIterationCondition(loop_nr, _opts.max_loops);
	old_residual = residual;

	// submit loop tuple
	_info.addTuple(IterationInformation::LoopTable);

	while (!stop_condition) {
		// update loop count
		++loop_nr;
		_info.replace(IterationInformation::LoopTable, "loop_nr", (int)loop_nr);

		BOOST_LOG_TRIVIAL(debug)
			<< "======================== #" << loop_nr << "/1 ==================";

//		renormalizeMatrixByTrace(spectral_matrix);

		if ((spectral_matrix.innerSize() > 10) || (spectral_matrix.outerSize() > 10)) {
			BOOST_LOG_TRIVIAL(trace)
					<< "Current spectral matrix is\n" << spectral_matrix;
		} else {
			BOOST_LOG_TRIVIAL(info)
					<< "Current spectral matrix is\n" << spectral_matrix;
		}

#ifdef MPI_FOUND
		if (world.size() == 1) {
#endif
			for (unsigned int dim = 0; dim < _data.cols(); ++dim) {
				Eigen::VectorXd solution;
				stop_condition &=
						!solve(_opts,
								spectral_matrix,
								_data.col(dim),
								pixel_matrix.col(dim),
								solution,
								dim,
								loop_nr);
				pixel_matrix.col(dim) = solution;
			}
#ifdef MPI_FOUND
		} else {
			stop_condition &=
					!solveMaster(
							world,
							_opts,
							spectral_matrix,
							_data,
							pixel_matrix,
							loop_nr);
		}
#endif

		if ((pixel_matrix.innerSize() > 10) || (pixel_matrix.outerSize() > 10)) {
			BOOST_LOG_TRIVIAL(trace)
					<< "Resulting pixel matrix is\n" << pixel_matrix;
		} else {
			BOOST_LOG_TRIVIAL(info)
					<< "Resulting pixel matrix is\n" << pixel_matrix;
		}

		// check criterion
		{
			residual = calculateResidual(_data, spectral_matrix, pixel_matrix);
			BOOST_LOG_TRIVIAL(info)
				<< "#" << loop_nr << " 1/2, residual is " << residual;
		}

		BOOST_LOG_TRIVIAL(debug)
			<< "======================== #" << loop_nr << "/2 ==================";

//		renormalizeMatrixByTrace(pixel_matrix);

		if ((pixel_matrix.innerSize() > 10) || (pixel_matrix.outerSize() > 10)) {
			BOOST_LOG_TRIVIAL(trace)
					<< "Current pixel matrix is\n" << pixel_matrix;
		} else {
			BOOST_LOG_TRIVIAL(info)
					<< "Current pixel matrix is\n" << pixel_matrix;
		}

		// must transpose in place, as spectral_matrix.transpose() is const
		spectral_matrix.transposeInPlace();
#ifdef MPI_FOUND
		if (world.size() == 1) {
# endif
			for (unsigned int dim = 0; dim < _data.rows(); ++dim) {
				Eigen::VectorXd solution;
				stop_condition &=
						solve(_opts,
								pixel_matrix.transpose(),
								_data.row(dim).transpose(),
								spectral_matrix.col(dim),
								solution,
								dim,
								loop_nr);
				spectral_matrix.col(dim) = solution;
			}
#ifdef MPI_FOUND
		} else {
			stop_condition &=
					solveMaster(
							world,
							_opts,
							pixel_matrix.transpose(),
							_data.transpose(),
							spectral_matrix,
							loop_nr);
		}
#endif
		spectral_matrix.transposeInPlace();

		if ((spectral_matrix.innerSize() > 10) || (spectral_matrix.outerSize() > 10)) {
			BOOST_LOG_TRIVIAL(trace)
					<< "Resulting spectral matrix is\n" << spectral_matrix;
		} else {
			BOOST_LOG_TRIVIAL(info)
					<< "Resulting spectral matrix is\n" << spectral_matrix;
		}

		// remove ambiguity
		renormalizeMatrixProduct(spectral_matrix, pixel_matrix);

		// check criterion
		{
			residual = calculateResidual(_data, spectral_matrix, pixel_matrix);
			BOOST_LOG_TRIVIAL(info)
				<< "#" << loop_nr << " 2/2, residual is " << residual;
			_info.replace(IterationInformation::LoopTable, "residual", residual);
			stop_condition =
					checkRelativeResidualCondition(old_residual, residual, _opts.residual_threshold)
					|| checkResidualCondition(residual, _opts.residual_threshold)
					|| checkIterationCondition(loop_nr, _opts.max_loops);
			old_residual = residual;
		}

		// submit loop tuple
		_info.addTuple(IterationInformation::LoopTable);
	}
	if (loop_nr > _opts.max_loops)
		BOOST_LOG_TRIVIAL(error)
			<< "Maximum number of loops " << _opts.max_loops
			<< " exceeded, stopping iteration.";
	else
		BOOST_LOG_TRIVIAL(info)
			<< "Loop iteration performed " << loop_nr
			<< " times.";
	_info.replace(IterationInformation::OverallTable, "loops", (int)loop_nr);
	_info.replace(IterationInformation::OverallTable, "residual", residual);
	_info.addTuple(IterationInformation::OverallTable);
	
#ifdef MPI_FOUND
	sendTerminate(world);
#endif

	/// output solution
	{
		using namespace MatrixIO;
		if (!_opts.solution_factor_one_file.string().empty()) {
			std::ofstream ost(_opts.solution_factor_one_file.string().c_str());
			if (ost.good())
				try {
					ost << spectral_matrix;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully write first solution factor to file.\n";
				}
			else {
				std::cerr << "Failed to open " << _opts.solution_factor_one_file.string() << std::endl;
				return 255;
			}
		} else {
			std::cout << "No first solution factor file name given." << std::endl;
		}

		if (!_opts.solution_factor_two_file.string().empty()) {
			std::ofstream ost(_opts.solution_factor_two_file.string().c_str());
			if (ost.good())
				try {
					ost << pixel_matrix;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully write second solution factor to file.\n";
				}
			else {
				std::cerr << "Failed to open " << _opts.solution_factor_two_file.string() << std::endl;
				return 255;
			}
		} else {
			std::cout << "No second solution factor file name given." << std::endl;
		}

		if (!_opts.solution_product_file.string().empty()) {
			std::ofstream ost(_opts.solution_product_file.string().c_str());
			if (ost.good())
				try {
					ost << spectral_matrix * pixel_matrix;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully write solution product to file.\n";
				}
			else {
				std::cerr << "Failed to open " << _opts.solution_product_file.string() << std::endl;
				return 255;
			}
		} else {
			std::cout << "No solution product file name given." << std::endl;
		}
	}
	BOOST_LOG_TRIVIAL(debug)
		<< "Resulting first factor transposed is\n" << spectral_matrix.transpose();
	BOOST_LOG_TRIVIAL(debug)
		<< "Resulting second factor is\n" << pixel_matrix;

	const Eigen::MatrixXd product_matrix = spectral_matrix * pixel_matrix;
	if ((_data.innerSize() <= 10) && (_data.outerSize() <= 10)) {
		BOOST_LOG_TRIVIAL(debug)
			<< "Data matrix was\n" << _data;
		BOOST_LOG_TRIVIAL(debug)
			<< "Product matrix is\n" << product_matrix;
		BOOST_LOG_TRIVIAL(info)
			<< "Difference matrix is\n" << _data - product_matrix;
	}
	BOOST_LOG_TRIVIAL(info)
		<< "Norm of difference is " << (_data - product_matrix).norm();

	return 0;
}

int main(int argc, char **argv)
{
	/// start MPI
#ifdef MPI_FOUND
	  mpi::environment env(argc, argv);
	  mpi::communicator world;
#endif

	int returnstatus = 0;
#ifdef MPI_FOUND
	if (world.rank() == 0) {
#endif
		/// starting timing
		boost::chrono::high_resolution_clock::time_point timing_start =
				boost::chrono::high_resolution_clock::now();

		/// some required parameters
		MatrixFactorizerOptions opts;
		if (returnstatus == 0)
			returnstatus = parseOptions(argc, argv, opts);

		/// parse the matrices
		Eigen::MatrixXd data;
		if (returnstatus == 0)
			returnstatus = parseDataFile(opts.data_file.string(), data);

		/// create Database
		IterationInformation info(opts, data.innerSize(), data.outerSize());

		/// perform factorization
		if (returnstatus != 0) {
#ifdef MPI_FOUND
			sendTerminate(world);
#endif
		} else
			returnstatus = MatrixFactorization(opts, data, info
#ifdef MPI_FOUND
				, world
#endif
				);

#ifdef MPI_FOUND
		// exchange return status
		mpi::broadcast(world, returnstatus, 0);
#endif

		/// finish timing
		boost::chrono::high_resolution_clock::time_point timing_end =
				boost::chrono::high_resolution_clock::now();
		BOOST_LOG_TRIVIAL(info) << "The operation took "
				<< boost::chrono::duration<double>(timing_end - timing_start)
				<< ".";
#ifdef MPI_FOUND
	} else {
		// enter in solve() loop
		worker(world);

		// exchange return status
		mpi::broadcast(world, returnstatus, 0);
	}

	// End MPI
	MPI_Finalize ();
#endif

	/// exit
	return returnstatus;
}
