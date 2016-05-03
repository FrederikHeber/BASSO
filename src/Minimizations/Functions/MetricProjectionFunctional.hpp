/*
 * MetricProjectionFunctional.hpp
 *
 *  Created on: May 20, 2015
 *      Author: heber
 */

#ifndef METRICPROJECTIONFUNCTIONAL_HPP_
#define METRICPROJECTIONFUNCTIONAL_HPP_

#include "BassoConfig.h"

#include <vector>

#include "Minimizations/types.hpp"

class Mapping;
class Norm;

/** Functor to calculate MetricProjectionFunctional functional/distance.
 *
 */
class MetricProjectionFunctional
{
public:
	/** Constructor for MetricProjectionFunctional.
	 *
	 * @param _dualnorm norm object of dual space
	 * @param _J_q duality mapping from dual space to space
	 * @param _dualpower power of this duality mapping
	 * \param _U search space spanned by column vectors
	 */
	MetricProjectionFunctional(
			const Norm &_dualnorm,
			const Mapping &_J_q,
			const double _dualpower,
			const std::vector<SpaceElement_ptr_t> &_U);

	/** Constructor for MetricProjectionFunctional.
	 *
	 * @param _problem inverse problem with refs to norm and duality mapping
	 * \param _U search space spanned by column vectors
	 */
	MetricProjectionFunctional(
			const InverseProblem_ptr_t &_problem,
			const std::vector<SpaceElement_ptr_t> &_U
			);

	~MetricProjectionFunctional() {}

	/** Implements MetricProjectionFunctional functional.
	 *
	 * \param _t coefficients to column vectors of search space
	 * \param _dualx current dual of solution
	 */
	double operator()(
			const std::vector<double> &_t,
			const SpaceElement_ptr_t &_dualx
			) const;

	/** Implements MetricProjectionFunctional functional.
	 *
	 * \param _t coefficients to column vectors of search space
	 * \param _dualx current dual of solution
	 */
	std::vector<double> gradient(
			const std::vector<double> &_t,
			const SpaceElement_ptr_t &_dualx
			) const;

private:

	/** Updates the dual iterate defined by the distance between
	 * \a _dualx and the element the expansion coefficients in
	 * \a _t with internal search space \a U
	 * @param _t expansion coefficients
	 * @param _dualx other vector
	 */
	void updateDualIterate(
			const std::vector<double> &_t,
			const SpaceElement_ptr_t &_dualx
			) const;

	/** Function for calculating norms of a vector of SpaceElements
	 * to allow storing them in constant member variable.
	 *
	 * @param _dualnorm norm object of dual space
	 * @param _U given vector of vectors
	 * @return vector of their norms
	 */
	static const std::vector<double> calculateNorms(
			const Norm &_dualnorm,
			const std::vector<SpaceElement_ptr_t> &_U
			);

private:
	//!> power type of the weight function of \a J_q
	const double dualpower;
	//!> lp Norm object
	const Norm &dualnorm;
	//!> Mapping object
	const Mapping &J_q;

	//!> search space for expanding with given coefficients
	const std::vector<SpaceElement_ptr_t> &U;
	//!> precomputed norms over all spanning vector of search space
	const std::vector<double> normsU;

	//! internal temporary variable
	mutable SpaceElement_ptr_t resx;

	//!> internal temporary vector for dual element
	mutable SpaceElement_ptr_t dual_resx;

	//!> internal vector with zero in all components
	const SpaceElement_ptr_t zeroVec;
};


#endif /* METRICPROJECTIONFUNCTIONAL_HPP_ */
