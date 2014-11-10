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

	/** Setter for useOptimalStepwidth.
	 *
	 * This is to have a definite place where useOptimalStepwidth is
	 * changed. Hence, it is const and cannot accidentally be changed in
	 * the code, but it can still be set after the instance has been created.
	 *
	 * @param _useOptimalStepwidth new value of useOptimalStepwidth
	 */
	void setuseOptimalStepwidth(const bool _useOptimalStepwidth);

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

private:
	/** Calculates tau for modulus of smoothness such that modulus
	 * over tau matches given \a _lambda
	 */
	double calculateMatchingTau(
			const SmoothnessModulus &_modul,
			const double _lambda) const;

public:
	/** Calculate optimal step width via a line search.
	 *
	 * i.e. along u we look for alpha to minimize residual
	 *
	 * @param _problem inverse problem with current solution
	 * @param _dualx dual element to current position
	 * @param _u direction for line search
	 * @return optimal step \f$ \alpha \f$ width to minimize
	 * 		\f$ || A (x - \alpha u) - y ||_Y \f$
	 */
	double calculateOptimalStepwidth(
			const InverseProblem_ptr_t &_problem,
			const SpaceElement_ptr_t &_dualx,
			const SpaceElement_ptr_t &_u,
			const double _alpha = 0) const;

public:
	//!> positive dampening constant for iteration
	const double C;
	//!> smoothness modulus object for dual Space X^*
	const SmoothnessModulus modul;
	//!> whether to use optimal step width calculation or a theoretical one
	const bool useOptimalStepwidth;
};


#endif /* LANDWEBERMINIMIZER_HPP_ */
