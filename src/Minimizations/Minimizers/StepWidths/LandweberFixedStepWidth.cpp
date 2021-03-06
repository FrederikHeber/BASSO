/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 * Copyright (C)  2016-2019 Frederik Heber. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*
 * LandweberFixedStepWidth.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "LandweberFixedStepWidth.hpp"

#include <boost/bind.hpp>
#include <boost/log/trivial.hpp>
#include <boost/math/tools/minima.hpp>
#include <limits>

#include "Log/Logging.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

LandweberFixedStepWidth::LandweberFixedStepWidth(
		const InverseProblem_ptr_t &_problem,
		const double _C) :
	problem(_problem),
	modul(_problem->x->getSpace()->getDualSpace()->getDualityMapping()->getPower()),
	C(_C)
{
	// some refs
	const Norm & NormX = *_problem->A->getSourceSpace()->getNorm();
	const Norm & DualNormX = *_problem->A->getSourceSpace()->getNorm();

	G = NormX.getPvalue() < 2. ?
			::pow(2., 2. - DualNormX.getPvalue()) :
			 DualNormX.getPvalue() - 1.;
	LOG(trace, "G is " << G);

	modulus_at_one = modul(1);
	ANorm = dynamic_cast<const LinearMapping &>(*_problem->A).Norm(); //::pow(2, 1.+ 1./val_NormY);
	LOG(trace, "ANorm " << ANorm);
}

const double LandweberFixedStepWidth::operator()(
		const SpaceElement_ptr_t &_dualx,
		const SpaceElement_ptr_t &_u,
		const SpaceElement_ptr_t &_solution,
		const SpaceElement_ptr_t &_residual,
		const double _residuum,
		const double _TolX,
		const double _alpha
		) const
{
	// some refs
	const Norm & NormX = *problem->x->getSpace()->getNorm();
	const Norm & DualNormX = *problem->x->getSpace()->getDualSpace()->getNorm();
	const Norm & NormY = *problem->y->getSpace()->getNorm();

	double alpha = _alpha;
	if (_dualx->isApproxToConstant(0, _TolX)) {
		const double q_p = ::pow(DualNormX.getPvalue(), NormX.getPvalue() - 1.);
		LOG(trace, "q_p " << q_p);
		const double A_p = ::pow(ANorm,NormX.getPvalue());
		LOG(trace, "A_p " << A_p);
		double R_p = 0.;
		if (NormX.getPvalue() == std::numeric_limits<double>::infinity()) {
			R_p = _residual->getMaxCoefficientAndIndex().first;
		} else {
			R_p = ::pow(_residuum, NormX.getPvalue());
		}
		double R_r = 0.;
		if (NormY.getPvalue() == std::numeric_limits<double>::infinity()) {
			R_r = _residual->getMaxCoefficientAndIndex().first;
		} else {
			R_r = ::pow(_residuum, NormY.getPvalue());
		}
		LOG(trace, "R_p " << R_p << ", R_r " << R_r);
		alpha = // C
				0.9 * (q_p / A_p) * (R_p/R_r);
	} else {
		const double xnorm = NormX(_solution);
		LOG(trace, "xnorm " << xnorm);
		const double two_q = ::pow(2., DualNormX.getPvalue());
		LOG(trace, "two_q " << two_q);
		alpha = C/(two_q * G * ANorm)
				* (_residuum / xnorm);
		LOG(trace, "initial lambda " << alpha);
		const double lambda = std::min(modulus_at_one, alpha);
		// find intermediate value in smoothness modulus to match lambda
		const double tau = calculateMatchingTau(lambda);
		// calculate step width
		const double x_p = ::pow( xnorm, NormX.getPvalue()-1.);
		LOG(trace, "x_p " << x_p);
		double R_r = 0.;
		if (NormY.getPvalue() == std::numeric_limits<double>::infinity()) {
			R_r = _residual->getMaxCoefficientAndIndex().first / _residuum;
		} else {
			R_r = ::pow(_residuum, NormY.getPvalue() - 1.);
		}
		LOG(trace, "R_r " << R_r);
		alpha = (tau/ANorm) * (x_p / R_r);
	}
	return alpha;
}

double LandweberFixedStepWidth::calculateMatchingTau(
		const double _lambda
		) const
{
	SmoothnessFunctional smoothness(modul, _lambda);
	double tau = .5;
	double minval = -1000.;
	double maxval = 1000.;
	int bits = 32;
	boost::uintmax_t maxiter = 100;
	std::pair<double, double> minpair =
			boost::math::tools::brent_find_minima(
					boost::bind(
							&SmoothnessFunctional::operator(),
							boost::cref(smoothness),
							_1),
					minval,
					maxval,
					bits,
					maxiter);
	tau = minpair.first;

	LOG(trace, "Matching tau from modulus of smoothness is " << tau);
	LOG(trace, "Counter-check: rho(tau)/tau = " << modul(tau)/tau << ", lambda = " << _lambda);

	return tau;
}
