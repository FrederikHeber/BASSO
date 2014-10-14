/*
 * SequentialSubspaceMinimizer.hpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#ifndef SEQUENTIALSUBSPACEMINIMIZER_HPP_
#define SEQUENTIALSUBSPACEMINIMIZER_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/GeneralMinimizer.hpp"

class Database;

/** This class implements the sequential subspace optimization by [Schöpfer,
 * Schuster,Louis, 2006].
 *
 */
class SequentialSubspaceMinimizer : public GeneralMinimizer
{
public:
	SequentialSubspaceMinimizer(
			const double _NormX,
			const double _NormY,
			const double _PowerX,
			const double _PowerY,
			const double _Delta,
			const unsigned int _maxiter,
			Database &_database,
			const unsigned int _outputsteps
			);
	~SequentialSubspaceMinimizer() {}

	/** Setter for N, the number of search directions.
	 *
	 * This is to have a definite place where N is changed. Hence,
	 * it is const and cannot accidentally be changed in the code, but
	 * it can still be set after the instance has been created.
	 *
	 * @param _N new value of N, N in [1, infty)
	 */
	void setN(const unsigned int _N);

	GeneralMinimizer::ReturnValues operator()(
			const Eigen::VectorXd &_x0,
			const Eigen::VectorXd &_dualx0,
			const Eigen::MatrixXd &_A,
			const Eigen::VectorXd &_y,
			const Eigen::VectorXd &_solution
			);

	/** Resets the iteration state of this minimizer in case
	 * the same object is to be used for another minimization with
	 * different problem matrix, right-hand side, ...
	 */
	void resetState()
	{ state.reset(); }

protected:
	/** This class encapsulates the state of the iteration, i.e. all
	 * variables that defined the current state of the iterate.
	 */
	struct IterationState
	{
		/** Default cstor for IterationState.
		 *
		 */
		IterationState() :
			isInitialized(false),
			index(-1)
		{}

		/** Resets the state to not initialized.
		 *
		 */
		void reset()
		{ isInitialized = false; }

		/** Initializer for the internal iteration state
		 *
		 * @param _dimension dimension of the search direction vectors
		 * @param _N number of search directions
		 */
		void set(
				const unsigned int _dimension,
				const unsigned int _N
				);

		/** Getter for the dimension of the search directions in \a U.
		 *
		 * @return size of search direction vectors,
		 * 		   0 - if state not initialized
		 */
		unsigned int getDimension() const
		{
			if (isInitialized)
				return U.innerSize();
			else
				return 0;
		}

		//!> flag to indicate whether inner state has been initialized
		bool isInitialized;
		//!> subspace matrix with search directions as column vectors
		Eigen::MatrixXd U;
		//!> offset of hyperplanes of search directions for projection
		Eigen::VectorXd alphas;
		//!> index of the last updated search direction
		unsigned int index;
	};

protected:
	// internal variables

//	//!> number of columns (M)
//	unsigned int NoCols;
//	//!> number of rows (N)
//	unsigned int NoRows;

	// constants

	//!> number of search directions
	const unsigned int N;

	//!> internal state of the iteration
	IterationState state;

	//!> counter for the small matrix vector products in subspace
	const OperationCounter<
		const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type,
		const Eigen::MatrixBase<Eigen::MatrixXd>&,
		const Eigen::MatrixBase<Eigen::VectorXd>&
		> MatrixVectorProduct_subspace;
	//!> counter for the small scalar products in subspace
	const OperationCounter<
		Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
		const Eigen::MatrixBase<Eigen::VectorXd>&,
		const Eigen::MatrixBase<Eigen::VectorXd>&
		> ScalarVectorProduct_subspace;
};


#endif /* SEQUENTIALSUBSPACEMINIMIZER_HPP_ */
