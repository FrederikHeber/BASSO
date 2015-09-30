/*
 * SpaceElementFactory.hpp
 *
 *  Created on: Jul 27, 2015
 *      Author: heber
 */

#ifndef SPACEELEMENTFACTORY_HPP_
#define SPACEELEMENTFACTORY_HPP_

#include "BassoConfig.h"

#include <string>

#include "Minimizations/types.hpp"

/** Either obtains matrix from a file or uses an internal one
 *
 */
class SpaceElementFactory
{
	/** Creates a SpaceElement by parsing a suitable vector from the
	 * file named \a _name.
	 *
	 * @param _SpaceRef reference to source space
	 * @param _name filename to parse vector from
	 * @return initialized element
	 */
	static SpaceElement_ptr_t create(
			const NormedSpace_weakptr_t _SpaceRef,
			const std::string &_name);

	/** Checks whether the string represents the path of a present file.
	 *
	 * @param _name name of file
	 * @return true file present, false - else
	 */
	static bool isPresentFile(const std::string &_name);
};



#endif /* SPACEELEMENTFACTORY_HPP_ */
