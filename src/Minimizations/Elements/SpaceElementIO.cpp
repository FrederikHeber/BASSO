/*
 * SpaceElementIO.cpp
 *
 *  Created on: Jul 27, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "SpaceElementIO.hpp"

#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

void SpaceElementWriter::output(std::ostream &_ost, const SpaceElement_ptr_t &_element)
{
	using namespace MatrixIO;
	_ost << RepresentationAdvocate::get(_element);
}

void SpaceElementReader::input(std::istream &_ist, SpaceElement_ptr_t &_element)
{
	using namespace MatrixIO;
	MatrixIO::dims_t dims = getDimensions(_ist);
	// assert its a row vector
	assert( _element->getSpace()->getDimension() == dims.first );
	assert( 1 == dims.second );
	Eigen::VectorXd vector;
	_ist >> vector;
	RepresentationAdvocate::set(_element, vector);
}
