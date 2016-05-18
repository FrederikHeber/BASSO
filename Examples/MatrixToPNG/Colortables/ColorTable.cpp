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
