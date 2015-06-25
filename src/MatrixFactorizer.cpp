
#include "BassoConfig.h"

#include <Eigen/Dense>
#include <fstream>
#include <string>

#ifdef MPI_FOUND
#include <mpi.h>
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
		const Eigen::MatrixXd &_rhs,
		Eigen::MatrixXd &_solution,
		const unsigned int _loop_nr
		)
{
	for (unsigned int dim = 0; dim < _rhs.cols(); ++dim) {
		Database_ptr_t mock_db(new Database_mock);
		Eigen::VectorXd projected_rhs_col(_rhs.col(dim));
		BOOST_LOG_TRIVIAL(debug)
			<< "------------------------ n=" << dim << " ------------------";
		// project right-hand side onto range of matrix
		BOOST_LOG_TRIVIAL(trace)
				<< "Initial y_" << dim << " is "
				<< projected_rhs_col.transpose();
		{
			// project y onto image of K
			if (!projectOntoImage(
					mock_db,
					_opts,
					_matrix,
					_rhs.col(dim),
					projected_rhs_col))
				return false;
		}
		BOOST_LOG_TRIVIAL(trace)
				<< "Projected y_" << dim << " is "
				<< projected_rhs_col.transpose();
		BOOST_LOG_TRIVIAL(trace)
				<< "Difference y_" << dim << "-y'_" << dim << " is "
				<< (_rhs.col(dim)-projected_rhs_col).transpose();
		BOOST_LOG_TRIVIAL(debug)
				<< "........................ n=" << dim << " ..................";

		// solve inverse problem for projected right-hand side
		Eigen::VectorXd solution_col(_solution.col(dim));
		if (!solveProblem(
				mock_db,
				_opts,
				_matrix,
				projected_rhs_col,
				_solution.col(dim),
				solution_col,
				_loop_nr >= 3))
			return false;
		_solution.col(dim) = solution_col;
		BOOST_LOG_TRIVIAL(trace)
			<< "Resulting vector is " << solution_col.transpose();
	}

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

int MatrixFactorization(
		const MatrixFactorizerOptions &_opts,
		const Eigen::MatrixXd &_data,
		IterationInformation &_info
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

		stop_condition &=
				solve(_opts,
						spectral_matrix,
						_data,
						pixel_matrix,
						loop_nr);

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
		stop_condition &=
				solve(_opts,
						pixel_matrix.transpose(),
						_data.transpose(),
						spectral_matrix,
						loop_nr);
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
	int my_rank = 0;

#ifdef MPI_FOUND
	int numprocs;
	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
#endif

	int returnstatus = 0;
	if (my_rank == 0) {
		int returnstatus = 0;
		/// starting timing
		boost::chrono::high_resolution_clock::time_point timing_start =
				boost::chrono::high_resolution_clock::now();

		/// some required parameters
		MatrixFactorizerOptions opts;
		returnstatus = parseOptions(argc, argv, opts);
		if (returnstatus != 0)
			return returnstatus;

		/// parse the matrices
		Eigen::MatrixXd data;
		returnstatus = parseDataFile(opts.data_file.string(), data);
		if (returnstatus != 0)
			return returnstatus;

		/// create Database
		IterationInformation info(opts, data.innerSize(), data.outerSize());

		/// perform factorization
		MatrixFactorization(opts, data, info);

		/// finish timing
		boost::chrono::high_resolution_clock::time_point timing_end =
				boost::chrono::high_resolution_clock::now();
		BOOST_LOG_TRIVIAL(info) << "The operation took "
				<< boost::chrono::duration<double>(timing_end - timing_start)
				<< ".";
	}
#ifdef MPI_FOUND
	// End MPI
	MPI_Finalize ();
#endif

	/// exit
	return returnstatus;
}
