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
			const DualityMappingsContainer &_container,
			const double _NormY,
			const double _PowerY,
			const double _Delta,
			const unsigned int _maxiter,
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
			const Eigen::VectorXd &_x0,
			const Eigen::VectorXd &_dualx0,
			const Eigen::MatrixXd &_A,
			const Eigen::VectorXd &_y,
			const Eigen::VectorXd &_solution
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
	 * @param _x current position
	 * @param _dualx dual element to current position
	 * @param _u direction for line search
	 * @param _A matrix
	 * @param _y right-hand side
	 * @param _alpha initial guess for step width
	 * @return optimal step \f$ \alpha \f$ width to minimize
	 * 		\f$ || A (x - \alpha u) - y ||_Y \f$
	 */
	double calculateOptimalStepwidth(
			const Eigen::VectorXd &_x,
			const Eigen::VectorXd &_dualx,
			const Eigen::VectorXd &_u,
			const Eigen::MatrixXd &_A,
			const Eigen::VectorXd &_y,
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
