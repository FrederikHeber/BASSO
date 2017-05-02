/*
 * NormFactory.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef NORMFACTORY_HPP_
#define NORMFACTORY_HPP_

#include "BassoConfig.h"

#include <map>
#include <string>
#include <vector>

#include <boost/any.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include "Minimizations/types.hpp"

class L1NormUnitTest;
class LInfinityNormUnitTest;
class LpNormUnitTest;

/** This factory instantiates the respective lp norm depending on the value
 * p.
 */
class NormFactory
{
	//!> grant L1Norm unit test access to "Space-less" creator function
	friend class L1NormUnitTest;
	//!> grant LInfinityNorm unit test access to "Space-less" creator function
	friend class LInfinityNormUnitTest;
	//!> grant LpNorm unit test access to "Space-less" creator function
	friend class LpNormUnitTest;
	//!> grant DualityMappingFactoryUnitTest unit test access to getMap()
	friend class DualityMappingFactoryUnitTest;

public:
	//!> typedef for the vector of arbitrary arguments.
	typedef std::vector<boost::any> args_t;

private:
	//!> typedef for bound functions that create a norm object
	typedef boost::function<
			Norm_ptr_t(
					const NormedSpace_weakptr_t,
					const args_t &)> NormCreator_t;

	/** typedef for the map to associated norm creator functions a with
	 * user comprehensible token.
	 */
	typedef std::map<std::string, NormCreator_t> TokenCreatorMap_t;

	/** Static function to return the filled map of type <-> creator
	 * associations.
	 *
	 * @return const ref to map instance
	 */
	static const TokenCreatorMap_t getMap();

	/** Creates an lp norm according to the given \a _p.
	 *
	 * @param _p p value of the lp norm
	 * @return Norm instance
	 */
	static Norm_ptr_t createLpInstance(const double _p);

public:

	/** Creates the norm of the desired type \a _token using the given
	 * arguments in \a _args.
	 *
	 * @param _token Token to identify type of norm
	 * @param _space space associated with this norm
	 * @param _args arguments for creating the norm
	 * @return
	 */
	static Norm_ptr_t create(
			const std::string &_token,
			const NormedSpace_weakptr_t _space,
			const args_t &_args);

	/** States whether the given type \a _token is valid, i.e. a creator
	 * for it is known.
	 *
	 * @param _token type to check
	 * @return true - creator present, false - else
	 */
	static bool isValidType(
			const std::string &_token
			);

private:
	/** Creates an lp norm according to the given \a _args.
	 *
	 * \sa LpNorm
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @param _args arguments containing exactly one double defining p.
	 * @return Norm instance
	 */
	static Norm_ptr_t createLpInstance(
			const NormedSpace_weakptr_t _ref,
			const args_t &_args);

	/** Creates the dual to an lp norm according to the given \a _args.
	 *
	 * \note The returned dual norm will have the conjugate value q with
	 * respect to p according to \f$ \frac 1 p + \frac 1 q = 1 \f$.
	 *
	 * \sa LpNorm
	 *
	 * @param _ref reference to the space this norm is associated with,
	 * 		  i.e. here the dual instance
	 * @param _args arguments containing exactly one double defining p.
	 * @return Norm instance
	 */
	static Norm_ptr_t createDualLpInstance(
			const NormedSpace_weakptr_t _ref,
			const args_t &_args);

	/** Creates a regularized l1 norm using relative shrinkage.
	 *
	 * \sa RelativeShrinkageL1Norm
	 *
	 * No arguments are needed. Lambda is set automatically depending on the
	 * argument.
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @param _args arguments containing exactly one double defining the
	 * 		  regularization parameter.
	 * @return Norm instance
	 */
	static Norm_ptr_t createRelativeShrinkageL1Instance(
			const NormedSpace_weakptr_t _ref,
			const args_t &_args);

	/** Creates a dual norm to regularized l1 norm using relative shrinkage.
	 *
	 * \sa DualRelativeShrinkageL1Norm
	 *
	 * No arguments are needed. Lambda is set automatically depending on the
	 * argument.
	 *
	 * @param _ref reference to the space this norm is associated with,
	 * 		  i.e. here the dual instance
	 * @param _args arguments containing exactly one double defining the
	 * 		  regularization parameter.
	 * @return Norm instance
	 */
	static Norm_ptr_t createDualRelativeShrinkageL1Instance(
			const NormedSpace_weakptr_t _ref,
			const args_t &_args);
	/** Creates a regularized l1 norm using soft thresholding.
	 *
	 * \sa SoftThresholdingL1Norm
	 *
	 * No arguments are needed. Lambda is set automatically depending on the
	 * argument.
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @param _args arguments containing exactly one double defining the
	 * 		  regularization parameter.
	 * @return Norm instance
	 */
	static Norm_ptr_t createSoftThresholdingL1Instance(
			const NormedSpace_weakptr_t _ref,
			const args_t &_args);

	/** Creates a dual norm to regularized l1 norm using soft thresholding.
	 *
	 * \sa DualSoftThresholdingL1Norm
	 *
	 * No arguments are needed. Lambda is set automatically depending on the
	 * argument.
	 *
	 * @param _ref reference to the space this norm is associated with,
	 * 		  i.e. here the dual instance
	 * @param _args arguments containing exactly one double defining the
	 * 		  regularization parameter.
	 * @return Norm instance
	 */
	static Norm_ptr_t createDualSoftThresholdingL1Instance(
			const NormedSpace_weakptr_t _ref,
			const args_t &_args);

	/** Creates an illegal norm.
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @return Norm instance
	 */
	static Norm_ptr_t createIllegalInstance(
			const NormedSpace_weakptr_t _ref,
			const args_t &_args);
};



#endif /* NORMFACTORY_HPP_ */
