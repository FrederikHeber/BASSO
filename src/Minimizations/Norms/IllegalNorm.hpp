/*
 * L1Norm.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef ILLEGALNORM_HPP_
#define ILLEGALNORM_HPP_

#include "BassoConfig.h"

#include <cmath>
#include <Eigen/Dense>

#include "Norm.hpp"

#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/NormExceptions.hpp"

/** This class implements an illegal norm that throws upon usage.
 *
 * This is to detect calls of norms in algorithms where the actual
 * norm is unknown and supposedly not required.
 *
 */
class IllegalNorm : public Norm
{
public:
//	/** Constructor of class L1Norm.
//	 *
//	 * @param _ref reference to the space this norm is associated with
//	 */
//	L1Norm(const boost::shared_ptr<NormedSpace> _ref) :
//		Norm(_ref)
//	{}

	const double operator()(const SpaceElement_ptr_t &_x) const
	{
		throw NormIllegalValue_exception()
				<< NormIllegalValue_name("IllegalNorm");
	}

	const double operator()(const Eigen::VectorXd &_x) const
	{
		throw NormIllegalValue_exception()
				<< NormIllegalValue_name("IllegalNorm");
	}
};


#endif /* ILLEGALNORM_HPP_ */
