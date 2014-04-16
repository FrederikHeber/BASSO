/*
 * LandweberMinimizer_calculateOptimalStepwidth.cpp
 *
 *  Created on: Apr 14, 2014
 *      Author: heber
 */

#include "LandweberMinimizer.hpp"

#define USE_APPROXIMATE_LINESEARCH

#include <Eigen/Dense>
#ifndef USE_APPROXIMATE_LINESEARCH
#include <levmar.h>
#else
#include <boost/bind.hpp>
#include <boost/math/tools/minima.hpp>
#endif
#include <limits>

#include "MinimizationExceptions.hpp"


struct MinimizationParameters
{
	MinimizationParameters(
		const Eigen::VectorXd &_x_n,
		const Eigen::VectorXd &_u_n,
		const Eigen::MatrixXd &_A,
		const Eigen::VectorXd &_y,
		const LandweberMinimizer &_landweber
			) :
		x_n(_x_n),
		u_n(_u_n),
		A(_A),
		y(_y),
		landweber(_landweber)
	{}

	const Eigen::VectorXd &x_n;
	const Eigen::VectorXd &u_n;
	const Eigen::MatrixXd &A;
	const Eigen::VectorXd &y;
	const LandweberMinimizer &landweber;
};

#ifndef USE_APPROXIMATE_LINESEARCH
/** Static function to calculate residual for given step width
 *
 */
static void
func_residual(double *p, double *hx, int m, int n, void *adata)
{
	struct MinimizationParameters *params =
			static_cast<MinimizationParameters *>(adata);
	// calculate new candidate position
	Eigen::VectorXd dual_solution =
			params->landweber.J_p(params->x_n, params->landweber.PowerX);
	dual_solution -= p[0] * params->u_n;
	Eigen::VectorXd x =
			params->landweber.J_q(dual_solution, params->landweber.DualPowerX);
	// calculate residual at candidate position (goes into hx[0])
	Eigen::VectorXd residual;
	hx[0] = params->landweber.calculateResidual(
			x,
			params->A,
			params->y,
			residual
			);
}
#endif

class function_residual
{
public:
	function_residual(
			const Eigen::VectorXd &_x,
			const Eigen::VectorXd &_u,
			const Eigen::MatrixXd &_A,
			const Eigen::VectorXd &_y,
			const LandweberMinimizer &_landweber
			) :
				x(_x),
				u(_u),
				A(_A),
				y(_y),
				landweber(_landweber)
	{}
	~function_residual() {}

	double operator()(double _arg) const
	{
		Eigen::VectorXd dual_solution =
				landweber.J_p(x, landweber.PowerX);
		dual_solution -= _arg * u;
		Eigen::VectorXd x =
				landweber.J_q(dual_solution, landweber.DualPowerX);
		// calculate residual at candidate position (goes into hx[0])
		Eigen::VectorXd residual;
		return landweber.calculateResidual( x, A, y, residual );
	}

private:
	const Eigen::VectorXd &x;
	const Eigen::VectorXd &u;
	const Eigen::MatrixXd &A;
	const Eigen::VectorXd &y;
	const LandweberMinimizer &landweber;
};

double LandweberMinimizer::calculateOptimalStepwidth(
		const Eigen::VectorXd &_x,
		const Eigen::VectorXd &_u,
		const Eigen::MatrixXd &_A,
		const Eigen::VectorXd &_y,
		const double _alpha) const
{
	double alpha = _alpha;
#ifndef USE_APPROXIMATE_LINESEARCH
	struct MinimizationParameters params(
			_x,	// x_n
			_u, // u_n
			_A, // A
			_y, // y
			*this // landweber
			);
	double info[LM_INFO_SZ];
	double x[] = { 0. };
	double opts[] = { TolX, TolX, TolX, TolX };
	const int m = 1;
	const int n = 1;
	double lb[] = { 0 };
	double ub[] = { std::numeric_limits<double>::max() };
	double *work = (double *)malloc( ( LM_DIF_WORKSZ(m, n) + m*m) * sizeof(double));
	double *covar = work+LM_DIF_WORKSZ(m, n);
	int ret = dlevmar_bc_dif(
			(*func_residual),
			&alpha,
			x,
			m,
			n,
			lb,
			ub,
			NULL,
			1000,
			opts, 	/* opts[4] */
			info,
			work,
			covar,
			&params
			);
	  if (ret == -1) {
		alpha = _alpha;
		throw MinimizationFunctionError_exception()
			<< MinimizationFunctionError_name(alpha);
	  }
	// free everything
	free(work);
#else
	function_residual res(
			_x,	// x_n
			_u, // u_n
			_A, // A
			_y, // y
			*this // landweber
			);
	double minval = -1000.;
	double maxval = 1000.;
	int bits = 32;
	boost::uintmax_t maxiter = 100;
	std::pair<double, double> minpair =
			boost::math::tools::brent_find_minima(
					boost::bind(
							&function_residual::operator(),
							boost::cref(res),
							_1),
					minval,
					maxval,
					bits,
					maxiter);
	alpha = minpair.first;
#endif

	return alpha;
}



