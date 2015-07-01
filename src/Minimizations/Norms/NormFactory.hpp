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
	 * @param _instance instance of NormFactory for creator bindings
	 * @return const ref to map instance
	 */
	static const TokenCreatorMap_t& getMap(
			const NormFactory &_instance);

	/** Creates an lp norm according to the given \a _p.
	 *
	 * @param _p p value of the lp norm
	 * @return Norm instance
	 */
	Norm_ptr_t createLpInstance(const double _p) const;

	/** Private default constructor for the class NormFactory.
	 *
	 * Instance should always be constructed via getInstance().
	 *
	 */
	NormFactory()
	{}

public:

	/** Returns the static instance to this singleton.
	 *
	 * @return static instance
	 */
	static const NormFactory& getInstance();

	/** Creates the norm of the desired type \a _token using the given
	 * arguments in \a _args.
	 *
	 * @param _token Token to identify type of norm
	 * @param _space space associated with this norm
	 * @param _args arguments for creating the norm
	 * @return
	 */
	Norm_ptr_t create(
			const std::string &_token,
			const NormedSpace_weakptr_t _space,
			const args_t &_args) const;

	/** States whether the given type \a _token is valid, i.e. a creator
	 * for it is known.
	 *
	 * @param _token type to check
	 * @return true - creator present, false - else
	 */
	bool isValidType(
			const std::string &_token
			) const;

private:
	/** Creates an lp norm according to the given \a _args.
	 *
	 * \sa LpNorm
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @param _args arguments containing exactly one double defining p.
	 * @return Norm instance
	 */
	Norm_ptr_t createLpInstance(
			const NormedSpace_weakptr_t _ref,
			const args_t &_args) const;

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
	Norm_ptr_t createDualLpInstance(
			const NormedSpace_weakptr_t _ref,
			const args_t &_args) const;

	/** Creates a regularized l1 norm.
	 *
	 * \sa RegularizedL1Norm
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @param _args arguments containing exactly one double defining the
	 * 		  regularization parameter.
	 * @return Norm instance
	 */
	Norm_ptr_t createRegularizedL1Instance(
			const NormedSpace_weakptr_t _ref,
			const args_t &_args) const;

	/** Creates a dual norm to regularized l1 norm.
	 *
	 * \sa DualRegularizedL1Norm
	 *
	 * @param _ref reference to the space this norm is associated with,
	 * 		  i.e. here the dual instance
	 * @param _args arguments containing exactly one double defining the
	 * 		  regularization parameter.
	 * @return Norm instance
	 */
	Norm_ptr_t createDualRegularizedL1Instance(
			const NormedSpace_weakptr_t _ref,
			const args_t &_args) const;

	/** Creates an illegal norm.
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @return Norm instance
	 */
	Norm_ptr_t createIllegalInstance(
			const NormedSpace_weakptr_t _ref,
			const args_t &_args) const;
};



#endif /* NORMFACTORY_HPP_ */
