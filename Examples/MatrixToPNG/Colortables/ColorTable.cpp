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
 * ColorTable.cpp
 *
 *  Created on: May 18, 2016
 *      Author: heber
 */


#include "BassoConfig.h"


#include "ColorTable.hpp"

ColorTable::ColorTable()
{
	// fill the table
	Eigen::MatrixXd redgreen(3,3);
	redgreen << 1, 0, 0,
			0, 0, 0,
			0, 1, 0;
	colortable.insert( std::make_pair ("redgreen", redgreen) );

	Eigen::MatrixXd redblue(3,3);
	redblue << 1, 0, 0,
			0, 0, 0,
			0, 0, 1;
	colortable.insert( std::make_pair ("redblue", redblue) );

	Eigen::MatrixXd bluegreenred(3,3);
	bluegreenred << 0, 0, 1,
			0, 1, 0,
			1, 0, 0;
	colortable.insert( std::make_pair ("bluegreenred", bluegreenred) );

	Eigen::MatrixXd blueblackred(3,3);
	blueblackred << 0, 0, 1,
			0, 0, 0,
			1, 0, 0;
	colortable.insert( std::make_pair ("blueblackred", blueblackred) );

	Eigen::MatrixXd blackwhite(2,3);
	blackwhite << 0, 0, 0,
			1, 1, 1;
	colortable.insert( std::make_pair ("blackwhite", blackwhite) );
}

const Eigen::MatrixXd& ColorTable::getColorTable( const std::string &_key)
{
	assert (colortable.count(_key) != 0);
	return colortable[_key];
}
