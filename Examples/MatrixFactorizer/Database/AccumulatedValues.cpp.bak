/*
 * AccumulatedValues.cpp
 *
 *  Created on: Nov 25, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "AccumulatedValues.hpp"

#include <iterator>

void AccumulatedValues::insert(
		const_iterator _first,
		const_iterator _last)
{
	for (const_iterator iter = _first; iter != _last; ++iter) {
		const std::string &keyname = iter->first;
		iterator insertiter = values.find(keyname);
		if (insertiter != values.end()) {
			// assure same types
			assert( iter->second.first == insertiter->second.first);
			// then insert all values
			insertiter->second.second.reserve(
					insertiter->second.second.size()+
					std::distance(iter->second.second.begin(), iter->second.second.end()));
			insertiter->second.second.insert(
					insertiter->second.second.end(),
					iter->second.second.begin(),
					iter->second.second.end());
		} else {
			std::pair<iterator, bool> inserter =
					values.insert( *iter );
			assert(inserter.second);
		}
	}
	values.insert(_first, _last);
}

size_t AccumulatedValues::getNumberOfValues() const
{
	if (values.size() == 0)
		return 0;
	else
		return values.begin()->second.second.size();
}
