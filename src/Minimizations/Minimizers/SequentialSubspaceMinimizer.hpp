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

#include "Minimizations/Minimizers/GeneralMinimizer.hpp"

class Database;

/** This class implements the sequential subspace optimization by [Schöpfer,
 * Schuster,Louis, 2006].
 *
 */
class SequentialSubspaceMinimizer : public GeneralMinimizer
{
public:
	SequentialSubspaceMinimizer(
			const InverseProblem_ptr_t &_inverseproblem,
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
			const InverseProblem_ptr_t &_problem,
			const SpaceElement_ptr_t &_startvalue,
			const SpaceElement_ptr_t &_dualstartvalue,
			const SpaceElement_ptr_t &_truesolution
			);

	/** Resets the iteration state of this minimizer in case
	 * the same object is to be used for another minimization with
	 * different problem matrix, right-hand side, ...
	 */
	void resetState()
	{ istate.reset(); }

protected:
	/** This class encapsulates the state of the iteration, i.e. all
	 * variables that define the current state of the iterate.
	 */
	struct IterationState
	{
		/** Default cstor for IterationState.
		 *
		 */
		IterationState();

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
				return U.outerSize();
			else
				return 0;
		}

		/** Method for updating the search space with a new direction.
		 *
		 * @param _newdir new search direction
		 * @param _iterate current iterate
		 * @param _alpha new offset
		 */
		void updateSearchSpace(
				const Eigen::VectorXd &_newdir,
				const Eigen::VectorXd &_iterate,
				const double _alpha);

		/** Const ref getter for \a U.
		 *
		 * @return const ref to U
		 */
		const Eigen::MatrixXd & getSearchSpace() const
		{ return U; }

		/** Const ref getter for \a alphas.
		 *
		 * @return const ref to alphas
		 */
		const Eigen::VectorXd & getAlphas() const
		{ return alphas; }

		/** Const ref getter for \a index.
		 *
		 * @return current index
		 */
		const unsigned int & getIndex() const {
			return index;
		}

	private:
		//!> subspace matrix with search directions as column vectors
		Eigen::MatrixXd U;
		//!> offset of hyperplanes of search directions for projection
		Eigen::VectorXd alphas;

		/** Function that simply advances the index in a round-robin fashion.
		 *
		 * @param _newdir new search direction
		 * @param _iterate current iterate
		 * @return index+1 mod N
		 */
		unsigned int advanceIndex(
				const Eigen::VectorXd &_newdir,
				const Eigen::VectorXd &_iterate) const;

	private:
		//!> flag to indicate whether inner state has been initialized
		bool isInitialized;

		//!> index of the last updated search direction
		unsigned int index;

		//!> bound function to allow various index update mechanisms
		typedef boost::function<
				unsigned int (
						const IterationState * const,
						const Eigen::VectorXd &,
						const Eigen::VectorXd &
						)> updater_t;
		updater_t updateIndex;
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
	IterationState istate;

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
