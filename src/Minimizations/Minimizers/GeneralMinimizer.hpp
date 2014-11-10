/*
 * GeneralMinimizer.hpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#ifndef GENERALMINIMIZER_HPP_
#define GENERALMINIMIZER_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>
#include <Eigen/Dense>

#include "MatrixIO/OperationCounter.hpp"

#include "Database/Table.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Functions/SmoothnessModulus.hpp"

class Database;

/** This class describes the interface to a general minimizer.
 *
 */
class GeneralMinimizer
{
public:
	GeneralMinimizer(
			const InverseProblem_ptr_t &_inverseproblem,
			const double _Delta,
			const unsigned int _maxiter,
			const unsigned int _maxinneriter,
			Database &_database,
			const unsigned int _outputsteps=0
			);

	virtual ~GeneralMinimizer() {}

	/** Internal structure for return values.
	 *
	 */
	struct ReturnValues
	{
		//!> solution vector
		SpaceElement_ptr_t m_solution;
		//!> residual vector
		SpaceElement_ptr_t m_residual;
		//!> remaining residuum, i.e. norm of residual
		double residuum;
		//!> number of outer iterations till solution
		int NumberOuterIterations;
	};

	/** Solve the inverse problem _A * x = _y for x with given
	 * tart value \a _x0, discretized operator \A _A and right-hand
	 * side \a _y.
	 *
	 * @param _x0 start value, may be zero vector
	 * @param _dualx0 dual element to initial \a _x0
	 * @param _A matrix as discretized operator
	 * @param _y right-hand side
	 * @param _solution additional true solution to calculate Bregman
	 * 			distance
	 * @return structure with solution and iteration information
	 */
	virtual ReturnValues operator()(
			const InverseProblem_ptr_t &_problem,
			const SpaceElement_ptr_t &_startvalue,
			const SpaceElement_ptr_t &_dualstartvalue,
			const SpaceElement_ptr_t &_truesolution
			) = 0;

	/** Resets the iteration state of this minimizer in case
	 * the same object is to be used for another minimization with
	 * different problem matrix, right-hand side, ...
	 */
	virtual void resetState() = 0;

	double calculateResidual(
			const InverseProblem_ptr_t &_problem,
			SpaceElement_ptr_t &_residual
			) const;

	/** Setter for MinLib via \a _name as string.
	 *
	 * @param _name name of library to use for minimization
	 */
	void setMinLib(const std::string &_name);

	/** Checks whether \a _name represents a valid name for a
	 * minimization library.
	 *
	 * @param _name name of library to check
	 * @return true - valid name, false - library name unknown
	 */
	bool isValidMinLibName(const std::string &_name);

	/** Checks whether walltime limit was surpassed.
	 *
	 * @param _time current used time (i.e. current time - start time)
	 * @return true - \a _time exceeds set walltime, false - may continue
	 */
	bool CheckWalltime(
			const boost::chrono::duration<double> &_time) const;

	/** Checks whether more iterations than a specified limit have
	 * been performed.
	 *
	 * \note This is overruled if \a MaxWalltime is given.
	 *
	 * @param _current_outeriterations number of iterations
	 * @return true - iteration count exceeds limit, false - may continue
	 */
	bool CheckIterations(
			const int _current_outeriterations) const;

	/** Checks whether the residuum has fallen beneath a specified
	 * limit.
	 *
	 * @param _residuum current residuum
	 * @return true - \a _residuum is less than \a TolY, false - may continue
	 */
	bool CheckResiduum(
			const double _residuum) const;

protected:
	/** Internal helper function for specific Minimizers to print debugging
	 *  solutions.
	 *
	 * @param _solution intermediate solution
	 * @param _A forward matrix, required for printing projected solution
	 * @param _NumberOuterIterations current iteration count
	 */
	void printIntermediateSolution(
			const Eigen::VectorXd &_solution,
			const Eigen::MatrixXd &_A,
			unsigned int _NumberOuterIterations
			) const;
public:
	//!> magnitude of noise
	const double Delta;
	//!> maximum time algorithm may spend on optimization
	boost::chrono::duration<double> MaxWalltime;
	//!> maximum number of iterations in outer loop
	const int MaxOuterIterations;
	//!> maximum number of iterations in inner loop
	const int MaxInnerIterations;
	//!> tolerance for objects in space X
	const double TolX;
	//!> tolerance for objects in space Y
	const double TolY;
	//!> tolerance for Fun
	const double TolFun;
	//!> output solution each .. steps, 0 means never
	unsigned int outputsteps;

	//!> enumeration of all available minimization libraries (for line search)
	enum MinimizationLibraries {
		gnuscientificlibrary=0,
		nonlinearoptimization=1,
		MAX_MinimizationLibraries
	};
	//!> stores which minimization library to use
	enum MinimizationLibraries MinLib;

	typedef std::map<std::string, enum MinimizationLibraries> MinLib_names_t;
	//!> name for each minimization library
	static MinLib_names_t MinLib_names;

protected:

	/** reference to an external database where we store infomation
	 * about the behavior of the iteration procedure.
	 */
	Database &database;

	static Table::Tuple_t preparePerIterationTuple(
			const double _val_NormX,
			const double _val_NormY,
			const unsigned int _N,
			const unsigned int _dim);

	static Table::Tuple_t prepareOverallTuple(
			const double _val_NormX,
			const double _val_NormY,
			const unsigned int _N,
			const unsigned int _dim);

	boost::function<
		const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type  (
				const Eigen::MatrixBase<Eigen::MatrixXd>&,
				const Eigen::MatrixBase<Eigen::VectorXd>&
				)> matrix_vector_fctor;
	const OperationCounter<
		const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type,
		const Eigen::MatrixBase<Eigen::MatrixXd>&,
		const Eigen::MatrixBase<Eigen::VectorXd>&
		> MatrixVectorProduct;

	boost::function<
		Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType (
				const Eigen::MatrixBase<Eigen::VectorXd>&,
				const Eigen::MatrixBase<Eigen::VectorXd>&)
				> scalar_vector_fctor;
	const OperationCounter<
		Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
		const Eigen::MatrixBase<Eigen::VectorXd>&,
		const Eigen::MatrixBase<Eigen::VectorXd>&
		> ScalarVectorProduct;
};


#endif /* GENERALMINIMIZER_HPP_ */
