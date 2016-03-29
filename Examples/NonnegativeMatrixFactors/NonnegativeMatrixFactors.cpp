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
			triplets.push_back( Eigen::Triplet<double>(i,j, _matrix.coeffRef(i,j)));
		}
	sparse_matrix.setFromTriplets( triplets.begin(), triplets.end() );

	return sparse_matrix;
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
	if (!opts.database_file.string().empty()) {
		database.setDatabaseFile(opts.database_file.string());
		Table& paramtable = database.addTable("parameters");
		Table::Tuple_t &paramtuple = paramtable.getTuple();
		paramtuple.insert( std::make_pair("truncation_dimension",
				(int)opts.truncation_dimension), Table::Parameter);

		paramtuple.insert( std::make_pair("cols", (int)matrix.cols()), Table::Parameter);
		paramtuple.insert( std::make_pair("rows", (int)matrix.rows()), Table::Parameter);
		paramtable.addTuple(paramtuple);

		size_t rowid = 1;
		if (database.isDatabaseFileGiven()) {
			// check for presence
			if (!database.isTuplePresentInTable(paramtable, paramtuple)) {
				BOOST_LOG_TRIVIAL(debug)
						<< "Parameter tuple not present, adding to table.";
				database.writeTable(paramtable);
			}
			// and store rowid
			rowid = database.getIdOfTuplePresentInTable(
					paramtable, paramtuple);
			// clear table such that present tuple is not stored again
			database.clearTable(paramtable.getName());
			BOOST_LOG_TRIVIAL(debug)
				<< "Setting parameter_key to " << rowid;
		}

		Table& datatable = database.addTable("data_overall");
		Table::Tuple_t &datatuple = datatable.getTuple();
		datatuple.insert( std::make_pair("parameters_fk", (int)rowid),
				Table::Parameter);
	}

	// starting timing
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	const unsigned int inner_dimension =
			std::min(opts.truncation_dimension,
					(unsigned int)std::min(matrix.rows(), matrix.cols()));
	if (inner_dimension != opts.truncation_dimension) {
		BOOST_LOG_TRIVIAL(warning)
				<< "Reduced truncation dimension to "
				<< inner_dimension << " because of maximum rank consideration.";
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
	BOOST_LOG_TRIVIAL(debug)
			<< "Singular values " << singularvalues;
	const Eigen::MatrixXd &Left_singularvectors = truncated_svd.getLeftSingularvectors();
	BOOST_LOG_TRIVIAL(debug)
			<< "Left eigenvectors " << Left_singularvectors;
	const Eigen::MatrixXd &Right_singularvectors = truncated_svd.getRightSingularvectors();
	BOOST_LOG_TRIVIAL(debug)
			<< "Right eigenvectors " << Right_singularvectors;
	// note that eigenvalues are given in increasing order
	BOOST_LOG_TRIVIAL(trace)
			<< "Matrix " << matrix;
	BOOST_LOG_TRIVIAL(trace)
			<< "U*S*V (truncated) " << Left_singularvectors * singularvalues.asDiagonal() * Right_singularvectors.transpose();
	// don't check for equality! It's only a truncated SVD!
//	assert( (USV - matrix).norm() < BASSOTOLERANCE);

	/// 2. Initialize W(:, 1) = sqrt(S(1, 1)) ∗ U(:, 1) and H(1, :) = sqrt(S(1, 1)) ∗ V(:, 1)'
	{
		// reversed order due to increasing eigenvalue sort order
		{
			const int index = 0;
			const double sv = sqrt(singularvalues(inner_dimension-1-index,0));
			W.col(index) = sv*Left_singularvectors.col(inner_dimension-1-index);
			H.row(index) = sv*Right_singularvectors.col(inner_dimension-1-index);
			BOOST_LOG_TRIVIAL(trace)
					<< "W.col(" << index << ")' " << W.col(index).transpose();
			BOOST_LOG_TRIVIAL(trace)
					<< "H.row(" << index << ") " << H.row(index);
		}

		/// for j = 2 : k
		for (int index = 1; index < (int)inner_dimension; ++index) {
			/// x =U(:,j); y =V(:,j);
			const Eigen::VectorXd x = Left_singularvectors.col(inner_dimension-1-index);
			const Eigen::VectorXd y = Right_singularvectors.col(inner_dimension-1-index);
			BOOST_LOG_TRIVIAL(trace)
					<< "#" << index << ": x " << x.transpose();
			BOOST_LOG_TRIVIAL(trace)
					<< "#" << index << ": y " << y.transpose();
			/// xp =pos(x); xn=neg(x); yp =pos(y); yn=neg(y);
			const Eigen::VectorXd xp = getPositiveSection(x);
			const Eigen::VectorXd xn = getNegativeSection(x);
			const Eigen::VectorXd yp = getPositiveSection(y);
			const Eigen::VectorXd yn = getNegativeSection(y);
			BOOST_LOG_TRIVIAL(trace)
					<< "#" << index << ": xp " << xp.transpose();
			BOOST_LOG_TRIVIAL(trace)
					<< "#" << index << ": xn " << xn.transpose();
			BOOST_LOG_TRIVIAL(trace)
					<< "#" << index << ": yp " << yp.transpose();
			BOOST_LOG_TRIVIAL(trace)
					<< "#" << index << ": yn " << yn.transpose();
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
			BOOST_LOG_TRIVIAL(trace)
					<< "W.col(" << index << ")' " << W.col(index).transpose();
			BOOST_LOG_TRIVIAL(trace)
					<< "H.row(" << index << ") " << H.row(index);
		}
		BOOST_LOG_TRIVIAL(trace)
				<< "matrix " << matrix;
		BOOST_LOG_TRIVIAL(trace)
				<< "A-W*H " << matrix-W*H;
		BOOST_LOG_TRIVIAL(debug)
				<< "|A - W*H|_2 " << (matrix-W*H).norm();
	}

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	const boost::chrono::duration<double> runtime =
			boost::chrono::duration<double>(timing_end - timing_start);
	BOOST_LOG_TRIVIAL(info) << "The operation took "
			<< runtime << ".";
	if (!opts.database_file.string().empty()) {
		Table& datatable = database.getTable("data_overall");
		Table::Tuple_t &datatuple = datatable.getTuple();
		datatuple.insert( std::make_pair("runtime", (double)runtime.count()), Table::Data);
		datatable.addTuple(datatuple);
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

