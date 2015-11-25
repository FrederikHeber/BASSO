/*
 * GeneralMinimizer.hpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#ifndef GENERALMINIMIZER_HPP_
#define GENERALMINIMIZER_HPP_

#include "BassoConfig.h"

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include "Database/Table.hpp"
#include "Minimizations/types.hpp"
#include "Minimizations/Minimizers/DatabaseManager.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion.hpp"

class BregmanDistance;
class CommandLineOptions;
class Database;
class IterationInformation;
class Mapping;
class QuickAccessReferences;

/** This class describes the interface to a general minimizer.
 *
 */
class GeneralMinimizer
{
public:
	GeneralMinimizer(
			const CommandLineOptions &_opts,
			const InverseProblem_ptr_t &_inverseproblem,
			Database &_database
			);

	virtual ~GeneralMinimizer();

	/** Internal structure containing the current search direction.
	 *
	 */
	struct SearchDirection {
		 SpaceElement_ptr_t Jw;
		 SpaceElement_ptr_t u;

		 /** Updates the search direction
		  *
		  * @param _refs quick access to Banach space objects
		  * @param _residual residual vector
		  */
		 void update(
		 	 	 const QuickAccessReferences &_refs,
		 	 	 const SpaceElement_ptr_t &_residual);
	} searchdir;

	/** Internal structure for return values.
	 *
	 */
	struct ReturnValues
	{
		//!> solution vector
		SpaceElement_ptr_t m_solution;
		//!> dual of solution vector
		SpaceElement_ptr_t m_dual_solution;
		//!> residual vector
		SpaceElement_ptr_t m_residual;
		//!> remaining residuum, i.e. norm of residual
		double residuum;
		//!> number of outer iterations till solution
		int NumberOuterIterations;

		//!> enumerating minimization states
		enum MinimizationStatus {
			notbegun=0,
			starting,
			started,
			error,
			finished
		};

		//!> status of the minimization
		enum MinimizationStatus status;

		/** Print the current state of return values.
		 *
		 * @param ynorm norm of y for relative residual
		 */
		void output(
				const double ynorm) const;
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
	 *
	 * \sa resetState_interal();
	 */
	void resetState();

	/** Calculates the residual for the given inverse \a _problem.
	 *
	 * @param _problem inverse problem containing matrix and right-hand side
	 * @param _residual residual vector on return
	 * @return magnitude of residual in the norm of its specified space
	 */
	double calculateResidual(
			const InverseProblem_ptr_t &_problem,
			SpaceElement_ptr_t &_residual
			) const;

	/** Helper function to calculate (and print) the Bregman distance.
	 *
	 * @param _Delta_p BregmanDistance object
	 * @param _solution current solution iterate
	 * @param _truesolution true solution
	 * @param _dual_solution dual of current solution iterate
	 * @return Bregman distance between \a _solution and \a _truesolution
	 */
	const double calculateBregmanDistance(
			const boost::shared_ptr<BregmanDistance> &_Delta_p,
			const SpaceElement_ptr_t &_solution,
			const SpaceElement_ptr_t &_truesolution,
			const SpaceElement_ptr_t &_dual_solution) const;

	/** Helper function to calculate (and print) the norm of the error.
	 *
	 * @param _solution current solution iterate
	 * @param _truesolution true solution
	 * @return norm of difference between \a _solution and \a _truesolution
	 */
	const double calculateError(
			const SpaceElement_ptr_t &_solution,
			const SpaceElement_ptr_t &_truesolution) const;


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

	/** Checks the stopping condition whether the iteration should stop
	 * or not.
	 *
	 * @param _time current used time (i.e. current time - start time)
	 * @param _current_outeriterations number of iterations
	 * @param _residuum current residuum
	 * @param _ynorm norm of the right-hand side
	 * @return true - stop, false - continue
	 */
	bool CheckStoppingCondition(
			const boost::chrono::duration<double> &_time,
			const int _current_outeriterations,
			const double _residuum,
			const double _ynorm) const;

	/** Sets the vector with additional parameters to make tuples unique in database.
	 *
	 * @param _tuple_params vector of strings in pairs of two ("key" and "value")
	 */
	void setAdditionalTupleParameters(
			const std::vector<std::string> &_tuple_params)
	{ dbcontainer.setAdditionalTupleParameters(_tuple_params); }

protected:
	/** Internal function called after GeneralMinimizer state has been
	 * resetted.
	 */
	virtual void resetState_interal()
	{}

	/** Resets the old Bregman distance if iteration is restarted or solver
	 * reused.
	 */
	void resetBregmanDistance() const
	{ OldBregmanDistance = 0.;	}

	/** Internal helper function for specific Minimizers to print debugging
	 *  solutions.
	 *
	 * @param _solution intermediate solution
	 * @param _A forward matrix, required for printing projected solution
	 * @param _NumberOuterIterations current iteration count
	 */
	void printIntermediateSolution(
			const SpaceElement_ptr_t &_solution,
			const Mapping &_A,
			unsigned int _NumberOuterIterations
			) const;

	/** Allows to add more parameters to the current parameter tuple before
	 * it is inserted into the sqlite database.
	 *
	 * @param _tuple tuple to insert parameters to
	 */
	virtual void addAdditionalParametersToTuple(
			Table::Tuple_t &_tuple,
			const bool _do_replace) const {}

protected:
	//!> IterationInformation needs access to DatabaseContainer
	friend class IterationInformation;

	//!> callback used by DatabaseManager needed in cstor
	DatabaseManager::addAdditionalParameters_t add_params_callback;

	//!> manager for all database access
	DatabaseManager dbcontainer;

public:
	//!> magnitude of noise
	const double Delta;
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

private:

	//!> contains the Bregman distance at the last iterate for monotonicity check
	mutable double OldBregmanDistance;

	//!> private l2 norm for calculating the error to the true solution
	const Norm_ptr_t l2norm;

protected:
	//!> stopping criteria
	const StoppingCriterion::ptr_t stopping_criteria;
};


#endif /* GENERALMINIMIZER_HPP_ */
