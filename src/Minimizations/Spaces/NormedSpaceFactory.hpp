/*
 * NormedSpaceFactory.hpp
 *
 *  Created on: Oct 27, 2014
 *      Author: heber
 */

#ifndef NORMEDSPACEFACTORY_HPP_
#define NORMEDSPACEFACTORY_HPP_

#include "BassoConfig.h"

#include <vector>

#include <string>
#include <boost/any.hpp>

#include "Minimizations/types.hpp"

class IllegalDualityMapping;
class NormFactory;

/** This is a factory for spaces, creating all necessary member instances
 * and returning the NormedSpace instance.
 *
 */
struct NormedSpaceFactory
{
	//!> typedef for the vector of arbitrary arguments.
	typedef std::vector<boost::any> args_t;

	/** Factory functor that creates the desired instance.
	 *
	 * @param _dimension dimension of the space
	 * @param _type desired type of space (and dual space)
	 * @param _args arguments required for instantiation
	 * @return NormedSpace instance according to parameters
	 */
	static NormedSpace_ptr_t create(
			const unsigned int _dimension,
			const std::string &_type,
			const args_t &_args);

private:
	//!> grant IllegalDualityMapping access to DummySpace
	friend class IllegalDualityMapping;
	//!> grant NormFactory access to DummySpace
	friend class NormFactory;

	/** This returns the ref to a placeholder instance to let entities
	 * such as DualityMappings receive a NormedSpace_ptr_t and to be
	 * able to call getDualSpace() properly.
	 *
	 * @return ref to placeholder instance
	 */
	static const NormedSpace_ptr_t getDummySpace();
};



#endif /* NORMEDSPACEFACTORY_HPP_ */
