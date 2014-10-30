/*
 * LandweberMinimizer.hpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#ifndef LANDWEBERMINIMIZER_HPP_
#define LANDWEBERMINIMIZER_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/StepWidths/DetermineStepWidthFactory.hpp"

class Database;

class LandweberMinimizer : public GeneralMinimizer
{
public:
	LandweberMinimizer(
			const InverseProblem_ptr_t &_inverseproblem,
			const double _Delta,
			const unsigned int _maxiter,
			const unsigned int _maxinneriter,
			Database &_database,
			const enum DetermineStepWidthFactory::stepwidth_enumeration _stepwidth_type,
			const unsigned int _outputsteps=0
			);

	~LandweberMinimizer() {}

	/** Setter for C.
	 *
	 * This is to have a definite place where C is changed. Hence,
	 * it is const and cannot accidentally be changed in the code, but
	 * it can still be set after the instance has been created.
	 *
	 * @param _C new value of C
	 */
	void setC(const double _C);

	GeneralMinimizer::ReturnValues operator()(
			const InverseProblem_ptr_t &_problem,
			const SpaceElement_ptr_t &_startvalue,
			const SpaceElement_ptr_t &_dualstartvalue,
			const SpaceElement_ptr_t &_truesolution
	);

	/** Resets the iteration state of this minimizer in case
	 * the same object is to be used for another minimization with
	 * different problem matrix, right-hand side, ...
	 *
	 * As the Landweber's inner state is empty, nothing needs to be done.
	 */
	void resetState() {}

public:
	//!> positive dampening constant for iteration
	const double C;
	//!> which step width procedure to use
	const enum DetermineStepWidthFactory::stepwidth_enumeration stepwidth_type;
};


#endif /* LANDWEBERMINIMIZER_HPP_ */
