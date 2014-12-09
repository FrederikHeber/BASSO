/*
 * VectorProjection_BregmanDistanceToLine.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef VECTORPROJECTION_BREGMANDISTANCETOLINE_HPP_
#define VECTORPROJECTION_BREGMANDISTANCETOLINE_HPP_

#include "BassoConfig.h"

#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/Functions/Minimizers/MinimizationFunctional.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/types.hpp"

/** This class implements a functional that measures the Bregman
 * Distance of a given vector to a line defined by another.
 */
class BregmanDistanceToLine :
		 public MinimizationFunctional<double>
{
	typedef typename MinimizationFunctional<double>::array_type array_type;
public:
	const BregmanDistance &distance;
	const Norm &dualnorm;
	const PowerTypeDualityMapping &J_q;
	//!> vector defining the line
	const SpaceElement_ptr_t &linevector;
	//!> vector who is to be projected onto linevector
	const SpaceElement_ptr_t &argvector;
	//!> power type of the distance
	const double powertype;
	//!> norm of line vector which is used multiply
	const double linevector_norm;

	BregmanDistanceToLine(
			const BregmanDistance &_distance,
			const Norm &_dualnorm,
			const PowerTypeDualityMapping &_J_q,
			const SpaceElement_ptr_t &_linevector,
			const SpaceElement_ptr_t &_tobeprojected,
			const double powertype
			);
	~BregmanDistanceToLine() {}

	double function(const double &_value) const;

	const double gradient(const double &_value) const;

	void convertInternalTypeToArrayType(
			const double &_t,
			array_type & _x
			) const;

	void convertArrayTypeToInternalType(
			const array_type & _x,
			double &_t
			) const;
};



#endif /* VECTORPROJECTION_BREGMANDISTANCETOLINE_HPP_ */
