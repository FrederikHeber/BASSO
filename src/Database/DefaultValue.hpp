/*
 * DefaultValue.hpp
 *
 *  Created on: Jun 13, 2014
 *      Author: heber
 */

#ifndef DEFAULTVALUE_HPP_
#define DEFAULTVALUE_HPP_

#include <string>

template <class T>
class DefaultValue {
public:
	DefaultValue() {}
	~DefaultValue() {};

	/** Getter for the contained value.
	 *
	 * @return returns the default value for the type T.
	 */
	T get() const
	{ return value; }

private:
	T value;
};

// specializations

template<>
inline
DefaultValue<int>::DefaultValue() :
	value(0)
{}

template<>
inline
DefaultValue<double>::DefaultValue() :
	value(0.)
{}

template<>
inline
DefaultValue<std::string>::DefaultValue() :
	value("-")
{}

#endif /* DEFAULTVALUE_HPP_ */
