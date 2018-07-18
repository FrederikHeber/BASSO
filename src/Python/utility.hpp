/* \file Utility - helpers/wrappers for pyBasso
 *
 * This file contains small wrapper functions to be used by pyBasso.
 *
 * utility.hpp
 *
 *  Created on: Jul 18, 2018
 *      Author: heber
 */

#ifndef PYTHON_UTILITY_HPP_
#define PYTHON_UTILITY_HPP_

#include <Eigen/Dense>

#include "Minimizations/types.hpp"

NormedSpace_ptr_t create_LpSpace(
		const int _dimension,
		const double _p,
		const double _power
		);

const Mapping_ptr_t create_LinearMapping(
		NormedSpace_ptr_t &_SourceSpaceRef,
		NormedSpace_ptr_t &_TargetSpaceRef,
		const Eigen::MatrixXd &_matrix,
		const bool _isAdjoint);

double SpaceElement_getitem(const SpaceElement_ptr_t &_element, const int _index);
void SpaceElement_setitem(SpaceElement_ptr_t &_element, const int _index, const double _value);

const bool SpaceElement_isZero(const SpaceElement_ptr_t &_element);
const bool SpaceElement_isApprox(
		const SpaceElement_ptr_t &_element,
		const SpaceElement_ptr_t &_other,
		const double _tolerance);

const Norm& NormedSpace_getNorm(const NormedSpace_ptr_t &_space);

const SpaceElement_ptr_t Mapping_operator(
		const Mapping_ptr_t &_map,
		const SpaceElement_ptr_t &_element);

const double Mapping_getTiming(const Mapping_ptr_t &_map);

#endif /* PYTHON_UTILITY_HPP_ */
