/*
 * GeneralMinimizer.hpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#ifndef GENERALMINIMIZER_HPP_
#define GENERALMINIMIZER_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "MatrixIO/OperationCounter.hpp"

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
	//!> maximum number of iterations in outer loop
	const int MaxOuterIterations;
	//!> tolerance for objects in space X
	const double TolX;
	//!> tolerance for objects in space Y
	const double TolY;
	//!> tolerance for Fun
	const double TolFun;
	//!> output solution each .. steps, 0 means never
	unsigned int outputsteps;

protected:

	/** reference to an external database where we store infomation
	 * about the behavior of the iteration procedure.
	 */
	Database &database;

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
