/*
 * SpaceElementFactory.cpp
 *
 *  Created on: Jul 27, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>

#include "MatrixIO/MatrixIO.hpp"
#include "MatrixIO/MatrixIOExceptions.hpp"
#include "Minimizations/Elements/SpaceElementFactory.hpp"
#include "Minimizations/Elements/SpaceElementIO.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

SpaceElement_ptr_t SpaceElementFactory::create(
		const NormedSpace_weakptr_t _SpaceRef,
		const std::string &_name
		)
{
	SpaceElement_ptr_t vector =
			NormedSpace_ptr_t(_SpaceRef)->createElement();

	using namespace MatrixIO;
	if (MatrixIO::isPresentFile(_name)) {
		std::ifstream ist(_name.c_str());
		if (ist.good()) {
			try {
				SpaceElementReader::input(ist, vector);
			} catch (MatrixIOStreamEnded_exception &e) {
				std::cerr << "Failed to fully parse vector from " << _name << std::endl;
			}
		} else {
			std::cerr << "Failed to open " << _name << std::endl;
		}
	} else {
		std::cerr << _name << " is not a valid name for SpaceElementFactory." << std::endl;
	}

	return vector;
}
