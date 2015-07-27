/*
 * SpaceElementIO.hpp
 *
 *  Created on: Jan 20, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_ELEMENTS_SPACEELEMENTIO_HPP_
#define MINIMIZATIONS_ELEMENTS_SPACEELEMENTIO_HPP_

#include "BassoConfig.h"

#include <iosfwd>

#include "Minimizations/types.hpp"

/** This is just an adaptor relating a SpaceElement to the MatrixIO.
 *
 * The reason is that we don't want MatrixIO to know about the specifics
 * of a SpaceElement and also we don't want to allow anyone to use
 * SpaceElement::getVectorRepresentation() which we need in order to
 * output the SpaceElement via MatrixIO. Hence, we wrap this access in
 * this static functor.
 */
struct SpaceElementWriter
{
	static void output(std::ostream &_ost, const SpaceElement_ptr_t &_element);
};

/** This is just an adaptor relating a SpaceElement to the MatrixIO.
 *
 * The reason is that we don't want MatrixIO to know about the specifics
 * of a SpaceElement and also we don't want to allow anyone to use
 * SpaceElement::getVectorRepresentation() which we need in order to
 * output the SpaceElement via MatrixIO. Hence, we wrap this access in
 * this static functor.
 */
struct SpaceElementReader
{
	static void input(std::istream &_ist, SpaceElement_ptr_t &_element);
};

#endif /* MINIMIZATIONS_ELEMENTS_SPACEELEMENTIO_HPP_ */
