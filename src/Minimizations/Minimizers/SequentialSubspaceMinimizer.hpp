/*
 * SequentialSubspaceMinimizer.hpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#ifndef SEQUENTIALSUBSPACEMINIMIZER_HPP_
#define SEQUENTIALSUBSPACEMINIMIZER_HPP_

#include "BassoConfig.h"

#include <set>
#include <vector>

#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/Searchspace/Searchspace.hpp"
#include "Minimizations/Minimizers/Searchspace/LastNSearchDirections.hpp"
#include "Minimizations/types.hpp"

class Database;

/** This class implements the sequential subspace optimization by [Schöpfer,
 * Schuster,Louis, 2006].
 *
 */
class SequentialSubspaceMinimizer : public GeneralMinimizer
{
public:
	SequentialSubspaceMinimizer(
			const CommandLineOptions &_opts,
			const InverseProblem_ptr_t &_inverseproblem,
			Database &_database
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
	 * For 1: Inexact line search will stop the iteration as soon as the Wolfe
	 * conditions are fulfilled.
	 * For 2: Inexact line search will stop the iteration as soon as
	 * a fixed count of iterations have been performed
	 *
	 * @param _inexactLinesearch: 0 - exact, 1 - inexact, 2 - fixed
	 */
	void setInexactLinesearch(const int _inexactLinesearch)
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

	/** Setter whether (bregman) angles between current and old search
	 * directions should be calculated or not.
	 *
	 * @param _flag true -  calculate angles, false - don't
	 */
	void setDoCalculateAngles(const bool _flag)
	{ DoCalculateAngles = _flag;}

	/** Resets the iteration state of this minimizer in case
	 * the same object is to be used for another minimization with
	 * different problem matrix, right-hand side, ...
	 */
	void resetState_interal()
	{ istate.reset(); }

	/** Sets the updateIndex function object of IterationState to
	 * the desired alternative update method.
 	 *
	 * @param _type desired variant of the updateIndex() algorithm
 	 */
	void setupdateIndexAlgorithm(
			const enum LastNSearchDirections::UpdateAlgorithmType _type);

protected:

	void addAdditionalParametersToTuple(
			Table::Tuple_t &_tuple,
			const bool _do_replace) const;

	const unsigned int calculateStepWidth(
			const QuickAccessReferences& refs,
			const SpaceElement_ptr_t& dual_solution,
			std::vector<double>& tmin,
			const std::vector<SpaceElement_ptr_t> &_searchspace,
			const std::vector<double> &_alphas
			) const;

	/** Checks whether we exceed \a initial_residuum by a certain
	 * factor and signals for stopping the iteration.
	 *
	 * \note Whether absolute or relative residuum is measured depends
	 * solely on the values given.
	 *
	 * @param current_residuum current residuum
	 * @param initial_residuum initial residuum
	 * @return true - stop iteration, false - continue
	 */
	bool isNonConverging(
			const double current_residuum,
			const double initial_residuum) const;

	void fillPerIterationTable(
			Table::Tuple_t& per_iteration_tuple,
			unsigned int tuple_counter);

	void fillAngleTable(
			Table::Tuple_t& angle_tuple,
			Table& data_angle_table,
			unsigned int tuple_counter);

	void updateAngleTable(
			const SpaceElement_ptr_t& newdir,
			Table::Tuple_t& angle_tuple) const;

	void updateSearchspace(
			const SpaceElement_ptr_t& _truesolution,
			const SpaceElement_ptr_t& newdir,
			const double alpha);

	void updateIterates(
			const QuickAccessReferences& refs,
			const std::vector<double> tmin,
			SpaceElement_ptr_t& _x,
			SpaceElement_ptr_t& dual_x) const;

	void updatePerIterationTuple(
			Table::Tuple_t& per_iteration_tuple,
			const double &ynorm,
			const double &stepwidth_norm,
			const int &inner_iterations,
			boost::shared_ptr<BregmanDistance> &Delta_p,
			const SpaceElement_ptr_t &_truesolution) const;

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
		void reset();

		/** Initializer for the internal iteration state
		 *
		 * @param _x0 initial iterate
		 * @param _dual_x0 initial dual of iterate
		 * @param _residual initial residual vector
		 * @param _residuum initial residuum
		 * @param _N number of search directions
		 * @param _orthogonalization_type orthogonalize new search direction?
		 */
		void set(
				const SpaceElement_ptr_t &_x0,
				const SpaceElement_ptr_t &_dual_x0,
				const SpaceElement_ptr_t &_residual,
				const double _residuum,
				const unsigned int _N,
				const LastNSearchDirections::OrthogonalizationType _orthogonalization_type
				);

		/** Getter for the dimension of the search directions in \a U.
		 *
		 * @return size of search direction vectors,
		 * 		   0 - if state not initialized
		 */
		unsigned int getDimension() const;

		/** Method for updating the search space with a new direction.
		 *
		 * @param _newdir new search direction
		 * @param _alpha new offset
		 * @param _dual_iterate current dual iterate
		 * @param _iterate current iterate
		 */
		void updateSearchSpace(
				const SpaceElement_ptr_t &_newdir,
				const double _alpha,
				const SpaceElement_ptr_t &_dual_iterate,
				const SpaceElement_ptr_t &_iterate
				);

		/** Const ref getter for representation of search space.
		 *
		 * @return const ref to U
		 */
		const std::vector<SpaceElement_ptr_t> & getSearchSpace() const;

		/** Const ref getter for representation of search space offsets.
		 *
		 * @return const ref to alphas
		 */
		const std::vector<double> & getAlphas() const;

		/** Getter for isInitialized, states whether set() needs to be
		 * called yet.
		 *
		 * @return true - set() was called, false - set() still needs to be called
		 */
		const bool getisInitialized() const;

		//!> typedef for a vector of angles
		typedef Searchspace::angles_t angles_t;

		/** Helper function to calculate the angles between each search
		 * direction in ::U and the given _newdir.
		 *
		 * @param _newdir new direction to compare to present ones
		 * @return vector of doubles, the angles
		 */
		const angles_t
		calculateAngles(
				const SpaceElement_ptr_t &_newdir) const;

		/** Helper function to calculate the angles between each search
		 * direction in ::U and the given _newdir using Bregman projections
		 * and distance.
		 *
		 * @param _newdir new direction to compare to present ones
		 * @return vector of doubles, the angles
		 */
		const angles_t
		calculateBregmanAngles(
				const SpaceElement_ptr_t &_newdir) const;

		//!> contains search space representation and update
		boost::shared_ptr<Searchspace> searchspace;

	private:
		/** Prevent copy cstor for IterationState.
		 *
		 */
		IterationState(const IterationState&);
		IterationState& operator=(const IterationState&);

	private:
		//!> flag to indicate whether inner state has been initialized
		bool isInitialized;
	};

protected:
	//!> number of search directions
	const unsigned int N;

	//!> internal state of the iteration
	IterationState istate;

	//!> bool whether to do an inexact line search with Wolfe conditions
	bool inexactLinesearch;

	//!> constant above which the step width must always lie
	double constant_positivity;
	//!> constant to scale the linear interpolation
	double constant_interpolation;

	//!> whether to calculate angles between current and old search dirs
	bool DoCalculateAngles;

	//!> whether to orthogonalize search directions w.r.t old ones
	const LastNSearchDirections::OrthogonalizationType OrthogonalizationType;

private:

	//!> internal temporary variable
	const SpaceElement_ptr_t dual_update;
};

inline
const bool SequentialSubspaceMinimizer::IterationState::getisInitialized() const
{ return isInitialized; }


#endif /* SEQUENTIALSUBSPACEMINIMIZER_HPP_ */
