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
 * printCounts.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef BASSO_PRINTCOUNTS_CPP_
#define BASSO_PRINTCOUNTS_CPP_

#include <Basso/printCounts.hpp>
#include "BassoConfig.h"


printCounts::printCounts(const std::vector<NormedSpace_ptr_t> &_spaces) :
	spaces(_spaces)
{}


#endif /* BASSO_PRINTCOUNTS_CPP_ */
