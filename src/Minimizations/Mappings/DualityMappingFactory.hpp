/*
 * DualityMappingFactory.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef DUALITYMAPPINGFACTORY_HPP_
#define DUALITYMAPPINGFACTORY_HPP_

#include "BassoConfig.h"

#include <map>
#include <string>
#include <vector>

#include <boost/any.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include "Minimizations/types.hpp"

struct DualityMappingFactory
{
public:
	//!> typedef for the vector of arbitrary arguments.
	typedef std::vector<boost::any> args_t;

private:
	//!> typedef for bound functions that create a norm object
	typedef boost::function<
			Mapping_ptr_t(
					const NormedSpace_weakptr_t,
					const args_t &)> MappingCreator_t;

	/** typedef for the map to associated norm creator functions a with
	 * user comprehensible token.
	 */
	typedef std::map<std::string, MappingCreator_t> TokenCreatorMap_t;

	/** Static function to return the filled map of type <-> creator
	 * associations.
	 *
	 * @param _instance instance of NormFactory for creator bindings
	 * @return const ref to map instance
	 */
	static const TokenCreatorMap_t& getMap(
			const DualityMappingFactory &_instance);

	/** Private default constructor for the class DualityMappingFactory.
	 *
	 * Instance should always be constructed via getInstance().
	 *
	 */
	DualityMappingFactory()
	{}

public:

	/** Returns the static instance to this singleton.
	 *
	 * @return static instance
	 */
	static const DualityMappingFactory& getInstance();

	/** Creates the norm of the desired type \a _token using the given
	 * arguments in \a _args.
	 *
	 * @param _token Token to identify type of norm
	 * @param _space space associated with this norm
	 * @param _args arguments for creating the norm
	 * @return
	 */
	Mapping_ptr_t create(
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
	/** Creates the duality mapping to an lp space.
	 *
	 * @param _NormedSpaceRef ref to normed space
	 * @param _args vector of arguments
	 * @return created instance
	 */
	Mapping_ptr_t createPowerTypeInstance(
			const NormedSpace_weakptr_t _NormedSpaceRef,
			const args_t &_args) const;

	/** Creates the duality mapping to an lp space.
	 *
	 * @param _NormedSpaceRef ref to normed space
	 * @param _args vector of arguments
	 * @return created instance
	 */
	Mapping_ptr_t createDualPowerTypeInstance(
			const NormedSpace_weakptr_t _NormedSpaceRef,
			const args_t &_args) const;

	/** Creates an illegal duality mapping.
	 *
	 * @param _NormedSpaceRef ref to normed space
	 * @param _args vector of arguments
	 * @return created instance
	 */
	Mapping_ptr_t createIllegalInstance(
			const NormedSpace_weakptr_t _NormedSpaceRef,
			const args_t &_args) const;

	/** Creates a relative shrinkage mapping.
	 *
	 * @param _NormedSpaceRef ref to normed space
	 * @param _args vector of arguments
	 * @return created instance
	 */
	Mapping_ptr_t createRelativeShrinkrageInstance(
			const NormedSpace_weakptr_t _NormedSpaceRef,
			const args_t &_args) const;
};



#endif /* DUALITYMAPPINGFACTORY_HPP_ */
