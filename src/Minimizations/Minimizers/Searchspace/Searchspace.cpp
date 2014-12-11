/*
 * Searchspace.cpp
 *
 *  Created on: Nov 17, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Searchspace.hpp"

#include <cmath>
#include <iterator>
#include <limits>

#include "Log/Logging.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Functions/VectorProjection.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

Searchspace::Searchspace(
		const NormedSpace_ptr_t &_SearchDirectionSpace_ptr,
		const unsigned int _N) :
	SearchDirectionSpace_ptr(_SearchDirectionSpace_ptr),
	projector(*_SearchDirectionSpace_ptr->getNorm(),
			dynamic_cast<const PowerTypeDualityMapping &>(
					*_SearchDirectionSpace_ptr->getDualityMapping()
					),
					_SearchDirectionSpace_ptr->getDualityMapping()->getPower()),
	U(_N),
	alphas(_N,0.)
{
	assert( U.size() == _N );
	std::generate(U.begin(), U.end(),
			boost::bind(&NormedSpace::createElement,
					boost::cref(_SearchDirectionSpace_ptr)));
	// DEBUG: check whether directions are initialized and zero
	const SearchDirections_t::const_iterator checkiter =
			std::find_if(U.begin(), U.end(),
					!boost::bind(&SpaceElement::isZero,
							_1,
							std::numeric_limits<double>::epsilon()*1e2));
	assert( checkiter == U.end() );
}

unsigned int
Searchspace::getDimension() const
{
	return getSearchSpace().size();
}

const Searchspace::angles_t
Searchspace::calculateBregmanAngles(
		const SpaceElement_ptr_t &_newdir) const
{
	// calculate all angles
	const unsigned int N = getDimension();
	angles_t angles(N, 0.);

	for (unsigned int l=0;l<N;++l) {
		if (U[l]->Norm() < std::numeric_limits<double>::epsilon())
			continue;
		// first: minimum, second: minimizer (equals length in LpNorm here)
		const std::pair<double, double> tmp =
				projector(
						U[l],
						_newdir,
						1e-4);
		const double projected_distance = tmp.second;
		const double original_distance = _newdir->Norm();
		if (fabs(original_distance) > std::numeric_limits<double>::epsilon()*1e2) {
			angles[l] = fabs(projected_distance / original_distance);
		} else {
			angles[l] = 0.;
		}

		BOOST_LOG_TRIVIAL(debug)
			<< "Bregman Angles #" << l << " is " << angles[l];
	}

	return angles;
}

const Searchspace::angles_t
Searchspace::calculateAngles(
		const SpaceElement_ptr_t &_newdir) const
{
	// calculate all angles
	const unsigned int N = getDimension();
	angles_t angles(N, 0.);

	for (unsigned int l=0;l<N;++l) {
		if (U[l]->Norm() < std::numeric_limits<double>::epsilon()) {
			continue;
		}
		// first: minimum, second: minimizer (equals length in LpNorm here)
		const double projected_distance = U[l] * _newdir / U[l]->Norm();
		const double original_distance = _newdir->Norm();
		if (fabs(original_distance) > std::numeric_limits<double>::epsilon()*1e2) {
			angles[l] = fabs(projected_distance / original_distance);
		} else {
			angles[l] = 0.;
		}
		BOOST_LOG_TRIVIAL(debug)
			<< "Angles #" << l << " is " << angles[l];
	}

	return angles;
}
