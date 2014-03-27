/*
 * DualityMapping.hpp
 *
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#ifndef DUALITYMAPPING_HPP_
#define DUALITYMAPPING_HPP_

#include "BassoConfig.h"

#include <cmath>
#include <Eigen/Dense>

#include "Math/Helpers.hpp"

/** This class contains a duality mapping instance from a specific
 * lp space to its dual.
 *
 * The template argument gives the value p of the lp norm.
 */
template <unsigned int p>
class DualityMapping
{
public:
	DualityMapping(const int _power) :
		m_power(_power)
	{}
	~DualityMapping() {}

	Eigen::VectorXd operator()(const Eigen::VectorXd &_x) const;

private:
	//!> contains the power of the weight (of power type) of the mapping
	const int m_power;
};

/** template specizialization for l_1 norm.
 *
 * \param _x vector
 * \return dual element corresponding to one element of the duality mapping for
 * 		x
 */
template <>
Eigen::VectorXd DualityMapping<1>::operator()(
		const Eigen::VectorXd &_x
		) const
{
	// J=norm(x,1)^(q-1)*sign(x);
	const double factor = ::pow(_x.lpNorm<1>(), (double)m_power-1.);
	return factor*Helpers::signum(_x);
}

/** template specizialization for l_infinity norm.
 *
 * \param _x vector
 * \return dual element corresponding to one element of the duality mapping for
 * 		x
 */
template <>
Eigen::VectorXd DualityMapping<Eigen::Infinity>::operator()(
		const Eigen::VectorXd &_x
		) const
{
	// [xNorm,k]=max(abs(x));
	unsigned int rowMax;
	unsigned int colMax;
	double factor = _x.array().abs().maxCoeff(&rowMax, &colMax);
	factor = ::pow(factor, (double)m_power-1.) * Helpers::sign(_x[rowMax]);
	// J=xNorm^(q-1)*sign(x(k,1))*circshift(eye(size(x)),[k-1 0]);
	Eigen::VectorXd temp = Eigen::VectorXd::Zero(_x.innerSize());
	temp[0] = 1.;
	const Eigen::VectorXd Jx = Helpers::circshift(temp, rowMax); // no -1 here, as index starts at 0 here, not 1
	return Jx * factor;
}

/** General function to calculate the duality mapping.
 *
 * \param _x vector
 * \return dual element corresponding to one element of the duality mapping for
 * 		x
 */
template <unsigned int p>
Eigen::VectorXd DualityMapping<p>::operator()(
		const Eigen::VectorXd &_x
		) const
{
	if (p == m_power) {
		// J=abs(x).^(p-1).*sign(x);
		Eigen::VectorXd Jx = _x.array().abs();
		for (int i=0;i<Jx.innerSize();++i)
			Jx[i] = ::pow(Jx[i], (double)p - 1.) * Helpers::sign(Jx[i]);
		return Jx;
	} else if ((int)p < m_power) {
		// J=norm(x,p)^(q-p)*abs(x).^(p-1).*sign(x);
		const double pnorm = ::pow(_x.lpNorm<p>(), (double)m_power-(double)p);
		Eigen::VectorXd Jx = _x.array().abs();
		for (int i=0;i<Jx.innerSize();++i)
			Jx[i] = pnorm * ::pow(Jx[i], (double)p - 1.) * Helpers::sign(Jx[i]);
		return Jx;
	} else {
		const double norm = _x.lpNorm<p>();
		if (norm < BASSOTOLERANCE) {
			// J=zeros(size(x,1),1);
			Eigen::VectorXd Jx = Eigen::VectorXd::Zero(_x.innerSize());
			return Jx;
		} else {
			// J=n^(q-p)*abs(x).^(p-1).*sign(x);
			const double exponent = (double)m_power-(double)p;
			const double pnorm = ::pow(norm, exponent);
			Eigen::VectorXd Jx = _x.array().abs();
			for (int i=0;i<Jx.innerSize();++i)
				Jx[i] = pnorm * ::pow(Jx[i], (double)p - 1.) * Helpers::sign(Jx[i]);
			return Jx;
		}
	}
}

#endif /* DUALITYMAPPING_HPP_ */
