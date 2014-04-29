/*
 * LandweberMinimizer_calculateMatchingTau.cpp
 *
 *  Created on: Apr 14, 2014
 *      Author: heber
 */

#include "LandweberMinimizer.hpp"

#include <boost/log/trivial.hpp>
#include <levmar.h>
#include <limits>

#include "MinimizationExceptions.hpp"

struct SmoothnessParameters
{
	SmoothnessParameters(
			const SmoothnessModulus &_modul,
			const double _lambda
			) :
		modul(_modul),
		lambda(_lambda)
	{}

	const SmoothnessModulus &modul;
	const double lambda;
};

/** Static function to calculate distance of smoothness modulus over tau
 * to given lambda with respect to tau.
 *
 */
static void
func_smoothness_over_tau(double *p, double *hx, int m, int n, void *adata)
{
	SmoothnessParameters *params = static_cast<SmoothnessParameters *>(adata);
	const double result = (params->modul)(p[0]);
	const double norm = result/p[0] - params->lambda;
	hx[0] = norm*norm;
}

double LandweberMinimizer::calculateMatchingTau(
		const SmoothnessModulus &_modul,
		const double _lambda
		) const
{
	SmoothnessParameters params(_modul, _lambda);
	double tau[] = { .5 };
	double info[LM_INFO_SZ];
	double x[] = { 0. };
	const int m = 1;
	const int n = 1;
	double lb[] = { std::numeric_limits<double>::epsilon() };
	double ub[] = { 1. };
	double *work = (double *)malloc( ( LM_DIF_WORKSZ(m, n) + m*m) * sizeof(double));
	double *covar = work+LM_DIF_WORKSZ(m, n);
	int ret = dlevmar_bc_dif(
			(*func_smoothness_over_tau),
			tau,
			x,
			m,
			n,
			lb,
			ub,
			NULL,
			1000,
			NULL, 	/* opts[4] */
			info,
			work,
			covar,
			&params
			);
//	  std::cout << "Minimization info: ";
//	  for(int i=0; i<LM_INFO_SZ; ++i)
//	    std::cout << info[i] << ", ";
//	  std::cout << std::endl;
	  if (ret == -1)
		throw MinimizationFunctionError_exception()
			<< MinimizationFunctionError_name(tau[0]);
	// free everything
	free(work);
	BOOST_LOG_TRIVIAL(trace)
		<< "Matching tau from modulus of smoothness is " << tau[0];
	BOOST_LOG_TRIVIAL(trace)
		<< "Counter-check: rho(tau)/tau = "
		<< _modul(tau[0])/tau[0] << ", lambda = " << _lambda;

	return tau[0];
}


