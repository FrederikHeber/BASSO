/*
 * L1Norm.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef L1NORM_HPP_
#define L1NORM_HPP_

#include "BassoConfig.h"

#include <cmath>
#include <Eigen/Dense>

#include "Norm.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"

/** This class implements the l_1 norm.
 *
 */
class L1Norm : public Norm
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
//		assert( NormedSpaceRef == x->getSpace() );
		return operator()(_x->getVectorRepresentation());
	}

	const double operator()(const Eigen::VectorXd &_x) const
	{
		double value = 0.;
		for (unsigned int i=0;i<_x.innerSize();++i)
			value += fabs(_x[i]);
		return value;
	}
};


#endif /* L1NORM_HPP_ */
