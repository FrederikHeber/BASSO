/*
 * BregmanProjectionFunctional.hpp
 *
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#ifndef BREGMANPROJECTIONFUNCTIONAL_HPP_
#define BREGMANPROJECTIONFUNCTIONAL_HPP_

#include "BassoConfig.h"

#include <vector>

#include "Minimizations/types.hpp"

class Mapping;
class Norm;

/** Functor to calculate BregmanProjectionFunctional functional/distance.
 *
 */
class BregmanProjectionFunctional
{
public:
	/** Constructor for BregmanProjectionFunctional.
	 *
	 * @param _dualnorm norm object of dual space
	 * @param _J_q duality mapping from dual space to space
	 * @param _dualpower power of this duality mapping
	 * \param _U search space spanned by column vectors
	 * \param _alpha offsets of affine subspace
	 */
	BregmanProjectionFunctional(
			const Norm &_dualnorm,
			const Mapping &_J_q,
			const double _dualpower,
			const std::vector<SpaceElement_ptr_t> &_U,
			const std::vector<double> &_alpha);

	/** Constructor for BregmanProjectionFunctional.
	 *
	 * @param _problem inverse problem with refs to norm and duality mapping
	 * \param _x current dual of solution
	 * \param _U search space spanned by column vectors
	 * \param _alpha offsets of affine subspace
	 */
	BregmanProjectionFunctional(
			const InverseProblem_ptr_t &_problem,
			const std::vector<SpaceElement_ptr_t> &_U,
			const std::vector<double> &_alpha);

	~BregmanProjectionFunctional() {}

	/** Implements BregmanProjectionFunctional functional.
	 *
	 * \param _t coefficients to column vectors of search space
	 * \param _dualx current dual of solution
	 */
	double operator()(
			const std::vector<double> &_t,
			const SpaceElement_ptr_t &dualx
			) const;

	/** Implements BregmanProjectionFunctional functional.
	 *
	 * \param _t coefficients to column vectors of search space
	 * \param _dualx current dual of solution
	 */
	std::vector<double> gradient(
			const std::vector<double> &_t,
			const SpaceElement_ptr_t &_dualx
			) const;

private:
	//!> power type of the weight function of \a J_q
	const double dualpower;
	//!> lp Norm object
	const Norm &dualnorm;
	//!> LpDualityMapping object
	const Mapping &J_q;

	//!> search space used for coefficients
	const std::vector<SpaceElement_ptr_t> &U;
	//!> offset to hyperplanes
	const std::vector<double> &alpha;
};


#endif /* BREGMANPROJECTIONFUNCTIONAL_HPP_ */
