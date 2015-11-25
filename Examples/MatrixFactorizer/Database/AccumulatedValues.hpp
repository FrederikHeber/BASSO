/*
 * AccumulatedValues.hpp
 *
 *  Created on: Nov 25, 2015
 *      Author: heber
 */

#ifndef ACCUMULATEDVALUES_HPP_
#define ACCUMULATEDVALUES_HPP_

#include "BassoConfig.h"

#include <map>
#include <string>

#include <boost/serialization/access.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>

#include "Database/Table.hpp"
#include "Database/types.hpp"

/** This represents a serializable version of the internal map to send
 * over network.
 */
struct AccumulatedValues
{
	//!> grant boost::serialization access to private members if any
	friend class boost::serialization::access;

	/** Serialization function for member variables.
	 *
	 * @param ar archive to store or load value in or from
	 * @param version version of to maintain compatibility
	 */
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & values;
	}

	typedef std::map<
				std::string,
				std::pair<
					Database_types::types_t,
					Table::any_values_t> > accumulatedValues_t;

	// expose some functionality of the internal map

	typedef accumulatedValues_t::iterator iterator;
	typedef accumulatedValues_t::const_iterator const_iterator;

	iterator find(const accumulatedValues_t::key_type &_key)
	{ return values.find(_key); }
	const_iterator find(const accumulatedValues_t::key_type &_key) const
	{ return values.find(_key); }

	accumulatedValues_t::mapped_type& operator[](
			const accumulatedValues_t::key_type &_key)
	{ return values[_key]; }

	void insert(
			const_iterator _first,
			const_iterator _last);

	accumulatedValues_t::size_type size() const
	{ return values.size(); }

	/** Returns the number of accumulated values for each key.
	 *
	 * @return single number, as all keys have same number of values.
	 */
	size_t getNumberOfValues() const;

	iterator end()
	{ return values.end(); }
	const_iterator end() const
	{ return values.end(); }

	iterator begin()
	{ return values.begin(); }
	const_iterator begin() const
	{ return values.begin(); }

	//!> typedef for internal values
	accumulatedValues_t values;
};


#endif /* ACCUMULATEDVALUES_HPP_ */
