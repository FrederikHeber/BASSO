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

#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/Specifics/LpDualityMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"

MetricProjectionFunctional::MetricProjectionFunctional(
		const Norm &_dualnorm,
		const PowerTypeDualityMapping &_J_q,
		const double _dualpower,
		const std::vector<SpaceElement_ptr_t> &_U
		) :
	dualpower(_dualpower),
	dualnorm(_dualnorm),
	J_q(_J_q),
	U(_U),
	normsU(calculateNorms(dualnorm,U))
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
	J_q(dynamic_cast<const PowerTypeDualityMapping&>(
			*_problem->SourceSpace->getDualSpace()->getDualityMapping())
			),
	U(_U),
	normsU(calculateNorms(dualnorm,U))
{
	assert( !U.empty() );
	assert( U.size() == normsU.size() );
}

SpaceElement_ptr_t MetricProjectionFunctional::calculateLinearCombination(
		const std::vector<double> &_t
		) const
{
	assert ( U.size() == _t.size() );
	SpaceElement_ptr_t resx = (*U.begin())->getSpace()->createElement();
	resx->setZero();
	*resx = std::inner_product(_t.begin(), _t.end(),U.begin(),resx);
	return resx;
}

SpaceElement_ptr_t MetricProjectionFunctional::calculateDistance(
		const std::vector<double> &_t,
		const SpaceElement_ptr_t &_dualx
		) const
{
	const SpaceElement_ptr_t resx = _dualx->getSpace()->createElement();
	*resx = _dualx;
	*resx -= calculateLinearCombination(_t);
	return resx;
}

double MetricProjectionFunctional::operator()(
		const std::vector<double> &_t,
		const SpaceElement_ptr_t &_dualx
		) const
{
	// x=x-U*t;
	const SpaceElement_ptr_t resx = calculateDistance(_t, _dualx);
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
	const SpaceElement_ptr_t dual_resx = J_q(calculateDistance(_t, _dualx));
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