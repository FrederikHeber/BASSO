/*
 * MetricProjectionFunctional.cpp
 *
 *  Created on: May 20, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include "MetricProjectionFunctional.hpp"

#include <boost/bind.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <numeric>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Norms/Norm.hpp"

MetricProjectionFunctional::MetricProjectionFunctional(
		const Norm &_dualnorm,
		const Mapping &_J_q,
		const double _dualpower,
		const std::vector<SpaceElement_ptr_t> &_U
		) :
	dualpower(_dualpower),
	dualnorm(_dualnorm),
	J_q(_J_q),
	U(_U),
	normsU(calculateNorms(dualnorm,U)),
	resx((*U.begin())->getSpace()->createElement()),
	dual_resx(resx->getSpace()->getDualSpace()->createElement()),
	zeroVec((*U.begin())->getSpace()->createElement())
{
	assert( !U.empty() );
	assert( U.size() == normsU.size() );
}

MetricProjectionFunctional::MetricProjectionFunctional(
		const InverseProblem_ptr_t &_problem,
		const std::vector<SpaceElement_ptr_t> &_U
		) :
	dualpower(_problem->SourceSpace->getDualSpace()->getDualityMapping()->getPower()),
	dualnorm(*_problem->SourceSpace->getDualSpace()->getNorm()),
	J_q(*_problem->SourceSpace->getDualSpace()->getDualityMapping()),
	U(_U),
	normsU(calculateNorms(dualnorm,U)),
	resx((*U.begin())->getSpace()->createElement()),
	dual_resx(resx->getSpace()->getDualSpace()->createElement()),
	zeroVec((*U.begin())->getSpace()->createElement())
{
	assert( !U.empty() );
	assert( U.size() == normsU.size() );
}

void MetricProjectionFunctional::updateDualIterate(
		const std::vector<double> &_t,
		const SpaceElement_ptr_t &_dualx
		) const
{
	resx->setZero();
	for (size_t dim = 0; dim < _t.size(); ++dim)
		resx->scaledAddition(-1.*_t[dim], U[dim]);
	*resx += _dualx;
}

double MetricProjectionFunctional::operator()(
		const std::vector<double> &_t,
		const SpaceElement_ptr_t &_dualx
		) const
{
	// x=x-U*t;
	updateDualIterate(_t, _dualx);
	// fval=1/q*norm(x,q)^q
	const double fval = 1./dualpower*::pow(dualnorm(resx), dualpower);
	return fval;
}

std::vector<double> MetricProjectionFunctional::gradient(
		const std::vector<double> &_t,
		const SpaceElement_ptr_t &_dualx
		) const
{
	// x=x-U*t;
	updateDualIterate(_t, _dualx);
	J_q(resx, dual_resx);
	std::vector<double> gval(_t.size(), 0.);
	for (size_t i=0;i<_t.size();++i)
		if (fabs(normsU[i] > BASSOTOLERANCE))
			gval[i] = -1.*(U[i] * dual_resx);
	return gval;
}

const std::vector<double> MetricProjectionFunctional::calculateNorms(
		const Norm &_dualnorm,
		const std::vector<SpaceElement_ptr_t> &_U
		)
{
	std::vector<double> norms(_U.size(), 0.);
	std::transform(
			_U.begin(), _U.end(),
			norms.begin(),
			boost::bind(&Norm::operator(),
					boost::cref(_dualnorm), _1));
	return norms;
}
