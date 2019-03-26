/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 * Copyright (C)  2016-2019 Frederik Heber. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

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
