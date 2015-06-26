/*
 * WorkResult.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_WORKRESULT_HPP_
#define MATRIXFACTORIZERBASE_WORKRESULT_HPP_

#include "BassoConfig.h"

#include <boost/serialization/access.hpp>

#include <Eigen/Dense>

/** Structure contains the result of a single work unit.
 *
 */
struct WorkResult
{
	/** Constructor for class WorkPackage for creating certain signaling
	 * work packages.
	 *
	 * \param _solution resulting solution of the inverse problem
	 * \param _solve_ok indicates whether solving was successful (true)
	 * \param _col column index unique to the work package
	 */
	WorkResult(
			const Eigen::VectorXd &_solution,
			const bool _solve_ok,
			const int _col
			) :
		solution(_solution),
		solve_ok(_solve_ok),
		col(_col)
	{}

	/** Default constructor for class WorkPackage.
	 *
	 */
	WorkResult() :
		solve_ok(true),
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
		ar & solution;
		ar & solve_ok;
		ar & col;
	}

	//!> solution of the inverse problem
	Eigen::VectorXd solution;
	//!> indicates whether a admissible solution has been obtained
	bool solve_ok;
	//!> column index unique to this work package
	int col;
};

#endif /* MATRIXFACTORIZERBASE_WORKRESULT_HPP_ */
