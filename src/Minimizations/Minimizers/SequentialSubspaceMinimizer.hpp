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
#include <set>
#include <vector>

#include "Minimizations/Functions/VectorProjection.hpp"
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

	/** Setter for whether the update algorithm is forced to be a random
	 * mapping.
	 *
	 * @param _enforceRandomMapping true - enforce update algorithm to be
	 *		  random mapping, false - no enforcement
	 */
	void setEnforceRandomMapping(const bool _enforceRandomMapping);

	/** Setter for whether an inexact or exact line search should be performed.
	 *
	 * Inexact line search will stop the iteration as soon as the Wolfe
	 * conditions are fulfilled.
	 *
	 * @param _inexactLinesearch true - do inexact line search, false - do not
	 */
	void setInexactLinesearch(const bool _inexactLinesearch)
	{ inexactLinesearch = _inexactLinesearch; }

	/** Setter for the constants in the three Wolfe conditions for the
	 * inexact line search.
	 *
	 * @param _wolfe_constants two values: positivity and stronger than linear
	 */
	void setWolfeConstants(const std::vector<double> &_wolfe_constants)
	{
		constant_positivity = _wolfe_constants[0];
		constant_interpolation = _wolfe_constants[1];
	}

	/** Resets the iteration state of this minimizer in case
	 * the same object is to be used for another minimization with
	 * different problem matrix, right-hand side, ...
	 */
	void resetState()
	{ istate.reset(); }

	//!> enumeration of available update index methods
	enum UpdateAlgorithmType {
		RoundRobin=0,
		MostParallel=1,
		MostOrthogonal=2,
		MAX_UpdateAlgorithmType
	};

	/** Sets the updateIndex function object of IterationState to
	 * the desired alternative update method.
 	 *
 	 * @param _problem inverse problem containing all instances
	 * @param _type desired variant of the updateIndex() algorithm
 	 */
	void setupdateIndexAlgorithm(
			const InverseProblem_ptr_t &_problem,
			const enum UpdateAlgorithmType _type);

protected:

	/** This class encapsulates the state of the iteration, i.e. all
	 * variables that define the current state of the iterate.
	 */
	struct IterationState : public ReturnValues
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
		 * @param _x0 initial iterate
		 * @param _residual initial residual vector
		 * @param _residuum initial residuum
		 * @param _N number of search directions
		 */
		void set(
				const SpaceElement_ptr_t &_x0,
				const SpaceElement_ptr_t &_residual,
				const double _residuum,
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
		 * @param _alpha new offset
		 */
		void updateSearchSpace(
				const SpaceElement_ptr_t &_newdir,
				const double _alpha);

		/** Const ref getter for \a U.
		 *
		 * @return const ref to U
		 */
		const Eigen::MatrixXd & getSearchSpace() const;

		/** Const ref getter for \a alphas.
		 *
		 * @return const ref to alphas
		 */
		const Eigen::VectorXd & getAlphas() const;

		/** Const ref getter for \a index.
		 *
		 * @return current index
		 */
		const unsigned int & getIndex() const {
			return index;
		}

		/** Getter for isInitialized, states whether set() needs to be
		 * called yet.
		 *
		 * @return true - set() was called, false - set() still needs to be called
		 */
		const bool getisInitialized() const;

		//!> typedef for a vector of angles
		typedef std::vector<double> angles_t;

		/** Helper function to calculate the angles between each search
		 * direction in ::U and the given _newdir.
		 *
		 * @param _Norm norm object to calculate norms
		 * @param _newdir new direction to compare to present ones
		 * @return vector of doubles, the angles
		 */
		const angles_t
		calculateAngles(
				const Norm &_Norm,
				const SpaceElement_ptr_t &_newdir) const;

		/** Helper function to calculate the angles between each search
		 * direction in ::U and the given _newdir using Bregman projections
		 * and distance.
		 *
		 * @param _Norm norm object to calculate norms
		 * @param _projector VectorProjection instance
		 * @param _newdir new direction to compare to present ones
		 * @return vector of doubles, the angles
		 */
		const angles_t
		calculateBregmanAngles(
				const Norm &_Norm,
				const VectorProjection &_projector,
				const SpaceElement_ptr_t &_newdir) const;

		//!> typedef for a set of indices
		typedef std::set<unsigned int> indexset_t;

	private:
		/** Prevent copy cstor for IterationState.
		 *
		 */
		IterationState(const IterationState&);
		IterationState& operator=(const IterationState&);

	private:
		//!> grant function access to updateIndex
		friend void SequentialSubspaceMinimizer::setupdateIndexAlgorithm(
				const InverseProblem_ptr_t &_problem,
				const enum UpdateAlgorithmType _type);
		//!> grant function access to setter for enforceRandomMapping
		friend void SequentialSubspaceMinimizer::setEnforceRandomMapping(
				const bool _enforceRandomMapping);

		//!> subspace matrix with search directions as column vectors
		Eigen::MatrixXd U;
		//!> offset of hyperplanes of search directions for projection
		Eigen::VectorXd alphas;

		/** Function that simply advances the index in a round-robin fashion.
		 *
		 * @param _newdir new search direction
		 * @return index+1 mod N
		 */
		unsigned int advanceIndex(
				const SpaceElement_ptr_t &_newdir) const;

		/** Setter for whether the update algorithm is forced to be a random
		 * mapping.
		 *
		 * @param _enforceRandomMapping true - enforce update algorithm to be
		 *		  random mapping, false - no enforcement
		 */
		void setEnforceRandomMapping(const bool _enforceRandomMapping)
		{  enforceRandomMapping = _enforceRandomMapping; }

		/** Function that compares new search direction with present ones
		 * in the state and selects the one that is most parallel,
		 * i.e. the one with maximum angle.
		 *
		 * @param _Norm norm object to calculate norms
		 * @param _projector VectorProjection instance
		 * @param _newdir new direction to compare to present ones
		 * @return index whose search direction is most parallel to \a newdir
		 */
		unsigned int
		updateIndexToMostParallel(
				const Norm &_Norm,
				const VectorProjection &_projector,
				const SpaceElement_ptr_t &_newdir) const;

		/** Function that compares new search direction with present ones
		 * in the state and selects the one that is most orthogonal,
		 * i.e. the one with minimum angle.
		 *
		 * @param _Norm norm object to calculate norms
		 * @param _projector VectorProjection instance
		 * @param _newdir new direction to compare to present ones
		 * @return index whose search direction is most orthogonal to \a newdir
		 */
		unsigned int
		updateIndexToMostOrthogonal(
				const Norm &_Norm,
				const VectorProjection &_projector,
				const SpaceElement_ptr_t &_newdir) const;

		/** Helper function for enforcing RandomMapping.
		 *
		 * This recreates the full index set.
		 *
		 * @param _indexset indexset to replenish
		 */
		void
		replenishIndexset(
				indexset_t &_indexset) const;

	private:
		//!> flag to indicate whether inner state has been initialized
		bool isInitialized;

		//!> index of the last updated search direction
		unsigned int index;

		//!> bound function to allow various index update mechanisms
		typedef boost::function<
				unsigned int (
						const IterationState * const,
						const SpaceElement_ptr_t &
						)> updater_t;
		updater_t updateIndex;

		//!> bool whether to enforce update index to be a random mapping or not
		bool enforceRandomMapping;

		//!> current set of indices remaining for enforcing random mapping
		mutable indexset_t current_indexset;
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

	//!> Vector Projection instance for calculating angles in Banach space
	VectorProjection projector;

	//!> bool whether to do an inexact line search with Wolfe conditions
	bool inexactLinesearch;

	//!> constant above which the step width must always lie
	double constant_positivity;
	//!> constant to scale the linear interpolation
	double constant_interpolation;
};

inline
const bool SequentialSubspaceMinimizer::IterationState::getisInitialized() const
{ return isInitialized; }


#endif /* SEQUENTIALSUBSPACEMINIMIZER_HPP_ */
