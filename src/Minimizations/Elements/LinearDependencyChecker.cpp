/*
 * LinearDependencyChecker.cpp
 *
 *  Created on: Feb 10, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "LinearDependencyChecker.hpp"

#include <Eigen/Dense>

#include "Log/Logging.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

bool LinearDependencyChecker::operator()(
		const vectors_t &_vectors)
{
	if (_vectors.empty()) {
		LOG(warning, "No vectors given whose linear dependency to check.");
		return true;
	}
	if (_vectors.size() == 1) {
		const bool result = _vectors[0]->isZero(BASSOTOLERANCE);
		BOOST_LOG_TRIVIAL(debug)
				<< "Checking Linear Dependency: Rank is "
				<< (result ? "0" : "1")
				<< " with " << _vectors.size() << " vectors.";
		return result;
	}
	Eigen::MatrixXd matrix(
			_vectors[0]->getSpace()->getDimension(),
			_vectors.size());
	for (size_t index = 0; index < _vectors.size(); ++index)
		matrix.col(index) =
				RepresentationAdvocate::get(_vectors[index]);
	const unsigned int rank = matrix.colPivHouseholderQr().rank();
	BOOST_LOG_TRIVIAL(debug)
			<< "Checking Linear Dependency: Rank is " << rank
			<< " with " << _vectors.size() << " vectors.";
	return rank != _vectors.size();
}


