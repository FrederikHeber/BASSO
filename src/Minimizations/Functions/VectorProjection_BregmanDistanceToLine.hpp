/*
 * VectorProjection_BregmanDistanceToLine.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef VECTORPROJECTION_BREGMANDISTANCETOLINE_HPP_
#define VECTORPROJECTION_BREGMANDISTANCETOLINE_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/Functions/MinimizationFunctional.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/Norm.hpp"

/** This class implements a functional that measures the Bregman
 * Distance of a given vector to a line defined by another.
 */
class BregmanDistanceToLine :
		 public MinimizationFunctional<double>
{
public:
	const BregmanDistance &distance;
	const Norm &dualnorm;
	const PowerTypeDualityMapping &J_q;
	//!> vector defining the line
	const Eigen::VectorXd &linevector;
	//!> vector who is to be projected onto linevector
	const Eigen::VectorXd &argvector;
	//!> power type of the distance
	const double powertype;

	BregmanDistanceToLine(
			const BregmanDistance &_distance,
			const Norm &_dualnorm,
			const PowerTypeDualityMapping &_J_q,
			const Eigen::VectorXd &_linevector,
			const Eigen::VectorXd &_tobeprojected,
			const double powertype
			);
	~BregmanDistanceToLine() {}

	double operator()(const double &_value) const;

	const double gradient(const double &_value) const;

	void convertToInternalType(
			double &_t,
			const gsl_vector * const x) const;

	void convertFromInternalType(
			const double &_t,
			gsl_vector * x) const;
};



#endif /* VECTORPROJECTION_BREGMANDISTANCETOLINE_HPP_ */
