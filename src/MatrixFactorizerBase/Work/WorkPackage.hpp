/*
 * WorkPackage.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_WORKPACKAGE_HPP_
#define MATRIXFACTORIZERBASE_WORKPACKAGE_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include <boost/serialization/access.hpp>

/** Structure contains all required to execute the work of a single
 * work package apart from GlobalData.
 */
struct WorkPackage
{
	/** Constructor for class WorkPackage.
	 *
	 * \note matrix is transfered as GlobalData
	 *
	 * @param _rhs right hand side
	 * @param _solution_startvalue starting value for solution
	 * @param _col column index (rhs and solution correspond to this column)
	 */
	WorkPackage(
			const Eigen::VectorXd &_rhs,
			const Eigen::VectorXd &_solution_startvalue,
			const int _col
			) :
		rhs(_rhs),
		solution_startvalue(_solution_startvalue),
		col(_col)
	{}

	/** Constructor for class WorkPackage for creating certain signaling
	 * work packages.
	 *
	 * \note i.e. column index of -1 indicates terminate.
	 *
	 * \param _col column index
	 */
	WorkPackage(
			const int _col) :
		col(_col)
	{}

	/** Default constructor for class WorkPackage.
	 *
	 */
	WorkPackage() :
		col(-1)
	{}

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
		ar & rhs;
		ar & solution_startvalue;
		ar & col;
	}

	//!> right hand side of the inverse problem
	Eigen::VectorXd rhs;
	//!> starting value for the sough-for solution
	Eigen::VectorXd solution_startvalue;
	//!> column index unique to the work package
	int col;
};





#endif /* MATRIXFACTORIZERBASE_WORKPACKAGE_HPP_ */
