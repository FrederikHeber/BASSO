/*
 * MetricProjectionFunctional.hpp
 *
 *  Created on: May 20, 2015
 *      Author: heber
 */

#ifndef METRICPROJECTIONFUNCTIONAL_HPP_
#define METRICPROJECTIONFUNCTIONAL_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>
#include <vector>

#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/types.hpp"

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
			const PowerTypeDualityMapping &_J_q,
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
	/** Helper function to calculate the element in the search space.
	 *
	 * @param _t expansion coefficients
	 * @return element defined by coefficients
	 */
	SpaceElement_ptr_t calculateLinearCombination(
			const std::vector<double> &_t
			) const;

	/** Calculates the distance vector between \a _dualx and the element
	 * defined by the expansion coefficients in \a _t with internal search
	 * space \a U
	 * @param _t expansion coefficients
	 * @param _dualx other vector
	 * @return distance vector
	 */
	SpaceElement_ptr_t calculateDistance(
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
	//!> LpDualityMapping object
	const PowerTypeDualityMapping &J_q;

	//!> search space for expanding with given coefficients
	const std::vector<SpaceElement_ptr_t> &U;
	//!> precomputed norms over all spanning vector of search space
	const std::vector<double> normsU;
};


#endif /* METRICPROJECTIONFUNCTIONAL_HPP_ */
