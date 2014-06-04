/*
 * LandweberMinimizer_calculateOptimalStepwidth.cpp
 *
 *  Created on: Apr 14, 2014
 *      Author: heber
 */

#include "LandweberMinimizer.hpp"

#include <Eigen/Dense>
#include <boost/bind.hpp>
#include <boost/math/tools/minima.hpp>
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

	return alpha;
}



