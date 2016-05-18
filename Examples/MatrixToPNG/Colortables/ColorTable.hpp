/*
 * ColorTable.hpp
 *
 *  Created on: May 18, 2016
 *      Author: heber
 */

#ifndef COLORTABLE_HPP_
#define COLORTABLE_HPP_

#include "BassoConfig.h"

#include <map>
#include <string>

#include <Eigen/Dense>

/** Color table serves as lookup for a number of predefine color table
 * matrices that can be used to colorize the data points.
 *
 * The color tables are stored as matrices of 3 columns representing
 * the RGB values. The rows designate the colors at regularly spaced
 * points between which we interpolate linearly.
 */
struct ColorTable
{
	/** Constructor for class ColorTable.
	 *
	 * Fills the internal table
	 *
	 */
	ColorTable();

	/** Checks for presence of \a _key in the table.
	 *
	 * @param _key key to check
	 * @return true - present, false - else
	 */	bool isKeyPresent(const std::string &_key)
	{ return colortable.count(_key); }

	 /** Getter for the color table to an existing \a _key.
	  *
	  * @return	Matrix representing the color table
	  */
	const Eigen::MatrixXd& getColorTable(const std::string &_key);

private:
	typedef std::map<std::string, Eigen::MatrixXd> colortable_t;

	//!> internal color table
	colortable_t colortable;
};


#endif /* COLORTABLE_HPP_ */
