/*
 * \file NonnegativeMatrixFactors.cpp
 *
 * This follows [Boutsidis, Gallopoulos, '08].
 *
 *  Created on: Nov 28, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <cmath>
#include <fstream>

#include <boost/chrono.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

#include "Database/SQLDatabase.hpp"
#include "Log/Logging.hpp"
#include "Options/NonnegativeMatrixFactorsOptions.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/Elements/SpaceElementIO.hpp"
#include "Solvers/SolverFactory/SolverFactory.hpp"
#include "Solvers/InverseProblemSolver.hpp"
#include "NonnegativeMatrixFactors/TruncatedSVD.hpp"

static Eigen::VectorXd getPositiveSection(const Eigen::VectorXd &_vector)
{
	Eigen::VectorXd returnvector = _vector.array().abs();
	returnvector += _vector;
	returnvector *= .5;
	return returnvector;
}

static Eigen::VectorXd getNegativeSection(const Eigen::VectorXd &_vector)
{
	Eigen::VectorXd returnvector = _vector.array().abs();
	returnvector -= _vector;
	returnvector *= .5;
	return returnvector;
}

Eigen::SparseMatrix<double> getSparseMatrix(const Eigen::MatrixXd& _matrix)
{
	// reserve space for full matrix
	Eigen::SparseMatrix<double> sparse_matrix(_matrix.rows(), _matrix.cols());
	Eigen::VectorXd nnz(_matrix.cols());
	for (int i=0;i<_matrix.cols();++i)
		nnz(i) = (double)_matrix.rows();
	sparse_matrix.reserve(nnz);

	// set all components
	std::vector< Eigen::Triplet<double> > triplets;
	triplets.reserve(_matrix.cols()*_matrix.rows());
	for (int i=0;i<_matrix.rows();++i)
		for (int j=0;j<_matrix.cols();++j) {
			if (fabs(_matrix.coeffRef(i,j)) > BASSOTOLERANCE)
				triplets.push_back( Eigen::Triplet<double>(i,j, _matrix.coeffRef(i,j)));
		}
	sparse_matrix.setFromTriplets( triplets.begin(), triplets.end() );

	return sparse_matrix;
}

size_t prepareDatabase(
		Database &_database,
		const NonnegativeMatrixFactorsOptions &_opts,
		const int _rows,
		const int _cols)
{
	_database.setDatabaseFile(_opts.database_file.string());
	Table& paramtable = _database.addTable("parameters");
	Table::Tuple_t &paramtuple = paramtable.getTuple();
	paramtuple.insert( std::make_pair("truncation_dimension",
			(int)_opts.truncation_dimension), Table::Parameter);
	paramtuple.insert( std::make_pair("cols", (int)_cols), Table::Parameter);
	paramtuple.insert( std::make_pair("rows", (int)_rows), Table::Parameter);
	paramtable.addTuple(paramtuple);

	size_t rowid = 1;
	if (_database.isDatabaseFileGiven()) {
		// check for presence
		if (!_database.isTuplePresentInTable(paramtable, paramtuple)) {
			LOG(debug, "Parameter tuple not present, adding to table.");
			_database.writeTable(paramtable);
		}
		// and store rowid
		rowid = _database.getIdOfTuplePresentInTable(
				paramtable, paramtuple);
		// clear table such that present tuple is not stored again
		_database.clearTable(paramtable.getName());
		LOG(debug, "Setting parameter_key to " << rowid);
	}

	return rowid;
}

void finalizeDatabase(
		Database &_database,
		const NonnegativeMatrixFactorsOptions &_opts)
{
	// write tables beforehand
	_database.writeAllTables();

	std::stringstream sql;
	sql << "CREATE VIEW IF NOT EXISTS overall AS SELECT * FROM parameters p INNER JOIN data_overall d ON p.rowid = d.parameters_fk";
	LOG(trace, "SQL: " << sql.str());
	if (!_database.executeSQLStatement(sql.str())) {
		LOG(error, "Failed to create view.");
	}
}

int main (int argc, char *argv[])
{
	// show program information
	showVersion(std::string(argv[0]));
	showCopyright();

	NonnegativeMatrixFactorsOptions opts;
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

	/// parse in source matrix
	Eigen::MatrixXd matrix;
	{
		using namespace MatrixIO;

		if (!MatrixIO::parse(opts.matrix.string(), "matrix factor", matrix))
			return 255;
	}

	// prepare database and parameter and data tables
	SQLDatabase database;
	size_t rowid = 0;
	if (!opts.database_file.string().empty())
		rowid = prepareDatabase(database, opts, matrix.cols(), matrix.rows());
	Table& datatable = database.addTable("data_overall");
	Table::Tuple_t &datatuple = datatable.getTuple();
	datatuple.insert( std::make_pair("parameters_fk", (int)rowid),
			Table::Parameter);

	// starting timing
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	const unsigned int inner_dimension =
			std::min(opts.truncation_dimension,
					(unsigned int)std::min(matrix.rows(), matrix.cols()));
	if (inner_dimension != opts.truncation_dimension) {
		LOG(warning, "Reduced truncation dimension to " << inner_dimension << " because of maximum rank consideration.");
	}

	/// apply method
	// create destination matrices
	Eigen::MatrixXd W(matrix.rows(), inner_dimension);
	Eigen::MatrixXd H(inner_dimension, matrix.cols());
	assert (W.cols() == H.rows());
	/// 1. Compute the largest k singular triplets of A
	TruncatedSVD::MatrixType_t sparse_matrix = getSparseMatrix(matrix);
	TruncatedSVD truncated_svd(sparse_matrix);
	truncated_svd(inner_dimension, TruncatedSVD::LargestMagnitude);
	const Eigen::MatrixXd &singularvalues = truncated_svd.getSingularvalues();
	LOG(debug, "Singular values of size " << singularvalues.rows() << "x" << singularvalues.cols() << " are\n" << singularvalues);
	const Eigen::MatrixXd &Left_singularvectors = truncated_svd.getLeftSingularvectors();
	LOG(debug, "Left eigenvectors of size " << Left_singularvectors.rows() << "x" << Left_singularvectors.cols() << " are\n" << Left_singularvectors);
	const Eigen::MatrixXd &Right_singularvectors = truncated_svd.getRightSingularvectors();
	LOG(debug, "Right eigenvectors of size " << Right_singularvectors.rows() << "x" << Right_singularvectors.cols() << " are\n" << Right_singularvectors);
	// note that eigenvalues are given in increasing order
	LOG(trace, "Matrix = \n" << matrix);
	LOG(debug, "| Matrix | = " << matrix.norm());
//	for (int i=0;i<opts.truncation_dimension;++i) {
//		Eigen::MatrixXd USV = Left_singularvectors.col(i) * singularvalues(i,0) * Right_singularvectors.col(i).transpose();
//		LOG(trace, "U*S*V (#" << i << ") = \n" << USV);
//	}
	Eigen::MatrixXd USV = Left_singularvectors * singularvalues.asDiagonal() * Right_singularvectors.transpose();
	const double Anorm = matrix.norm();
	const double USVnorm = USV.norm();
	const double SVDDifferencenorm = (matrix - USV).norm();
	LOG(trace, "U*S*V (truncated) = \n" << USV);
	LOG(debug, "|U*S*V (truncated)| = " << USVnorm);
	LOG(debug, "|U*S*V (truncated) - matrix| = " << SVDDifferencenorm);
	// don't check for equality! It's only a truncated SVD!
	// hence, we only check that the residual has gotten smaller
	assert( SVDDifferencenorm < Anorm);

	/// 2. Initialize W(:, 1) = sqrt(S(1, 1)) ∗ U(:, 1) and H(1, :) = sqrt(S(1, 1)) ∗ V(:, 1)'
	double residual = 0.;
	{
		// reversed order due to increasing eigenvalue sort order
		{
			const int index = 0;
			const double sv = sqrt(singularvalues(inner_dimension-1-index,0));
			W.col(index) = sv*Left_singularvectors.col(inner_dimension-1-index);
			H.row(index) = sv*Right_singularvectors.col(inner_dimension-1-index);
			LOG(trace, "W.col(" << index << ")' " << W.col(index).transpose());
			LOG(trace, "H.row(" << index << ") " << H.row(index));
		}

		/// for j = 2 : k
		for (int index = 1; index < (int)inner_dimension; ++index) {
			/// x =U(:,j); y =V(:,j);
			const Eigen::VectorXd x = Left_singularvectors.col(inner_dimension-1-index);
			const Eigen::VectorXd y = Right_singularvectors.col(inner_dimension-1-index);
			LOG(trace, "#" << index << ": x " << x.transpose());
			LOG(trace, "#" << index << ": y " << y.transpose());
			/// xp =pos(x); xn=neg(x); yp =pos(y); yn=neg(y);
			const Eigen::VectorXd xp = getPositiveSection(x);
			const Eigen::VectorXd xn = getNegativeSection(x);
			const Eigen::VectorXd yp = getPositiveSection(y);
			const Eigen::VectorXd yn = getNegativeSection(y);
			LOG(trace, "#" << index << ": xp " << xp.transpose());
			LOG(trace, "#" << index << ": xn " << xn.transpose());
			LOG(trace, "#" << index << ": yp " << yp.transpose());
			LOG(trace, "#" << index << ": yn " << yn.transpose());
			/// xpnrm=norm(xp); ypnrm=norm(yp);mp =xpnrm∗ ypnrm
			const double xp_norm = xp.norm();
			const double yp_norm = yp.norm();
			const double mp = xp_norm * yp_norm;
			/// xnnrm=norm(xn); ynnrm=norm(yn);mn=xnnrm∗ ynnrm;
			const double xn_norm = xn.norm();
			const double yn_norm = yn.norm();
			const double mn = xn_norm * yn_norm;
			assert( mp >= 0. );
			assert( mn >= 0. );
			/// if mp>mn,
			if ((mp > BASSOTOLERANCE) || (mn > BASSOTOLERANCE)) {
				const double sv = sqrt(singularvalues(inner_dimension-1-index,0));
				if (mp > mn) {
					/// u=xp/xpnrm; v =yp/ypnrm; sigma =mp;
					/// W(:,j)=sqrt(S(j, j) ∗ sigma) ∗ u and H(j, :)=sqrt(S(j, j) ∗ sigma) ∗ v?;
					const double smp = sqrt(mp);
					W.col(index) = sv*(smp/xp_norm)*xp;
					H.row(index) = sv*(smp/yp_norm)*yp;
					/// W(:,j)=sqrt(S(j, j) ∗ sigma) ∗ u and H(j, :)=sqrt(S(j, j) ∗ sigma) ∗ v?;
				} else {
					/// u=xn/xnnrm; v =yn/ynnrm; sigma =mn; end
					/// W(:,j)=sqrt(S(j, j) ∗ sigma) ∗ u and H(j, :)=sqrt(S(j, j) ∗ sigma) ∗ v?;
					const double smn = sqrt(mn);
					W.col(index) = sv*(smn/xn_norm)*xn;
					H.row(index) = sv*(smn/yn_norm)*yn;
				}
			} else {
				W.col(index).setZero();
				H.row(index).setZero();
			}
			LOG(trace, "W.col(" << index << ")' " << W.col(index).transpose());
			LOG(trace, "H.row(" << index << ") " << H.row(index));
		}
		LOG(trace, "matrix " << matrix);
		LOG(trace, "A-W*H " << matrix-W*H);
		residual = (matrix-W*H).norm();
		LOG(debug, "|A - W*H|_2 " << residual);
	}

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	const boost::chrono::duration<double> runtime =
			boost::chrono::duration<double>(timing_end - timing_start);
	LOG(info, "The operation took " << runtime << ".");

	// write runtime to database
	if (!opts.database_file.string().empty()) {
		Table& datatable = database.getTable("data_overall");
		Table::Tuple_t &datatuple = datatable.getTuple();
		datatuple.insert( std::make_pair("runtime", (double)runtime.count()), Table::Data);
		datatuple.insert( std::make_pair("Anorm", Anorm), Table::Data);
		datatuple.insert( std::make_pair("USVnorm", USVnorm), Table::Data);
		datatuple.insert( std::make_pair("SVDDifferencenorm", SVDDifferencenorm), Table::Data);
		datatuple.insert( std::make_pair("residual", residual), Table::Data);
		datatable.addTuple(datatuple);

		finalizeDatabase(database, opts);
	}

	/// store destination matrix factors
	if (!MatrixIO::store(
			opts.destination_first_factor.string(),
			"first matrix factor",
			W))
		return 255;
	if (!MatrixIO::store(
			opts.destination_second_factor.string(),
			"second matrix factor",
			H))
		return 255;

	// exit
	return 0;
}

