/*
 * types.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef TYPES_HPP_
#define TYPES_HPP_

#include "BassoConfig.h"

#include <boost/shared_ptr.hpp>

class InverseProblem;
class Mapping;
class Norm;
class NormedSpace;
class PowerTypeDualityMapping;
class SpaceElement;

typedef boost::shared_ptr<InverseProblem> InverseProblem_ptr_t;
typedef boost::shared_ptr<Mapping> Mapping_ptr_t;
typedef boost::shared_ptr<Norm> Norm_ptr_t;
typedef boost::shared_ptr<NormedSpace> NormedSpace_ptr_t;
typedef boost::shared_ptr<PowerTypeDualityMapping> PowerTypeDualityMapping_ptr_t;
typedef boost::shared_ptr<SpaceElement> SpaceElement_ptr_t;



#endif /* TYPES_HPP_ */