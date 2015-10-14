/*
 * StoppingArguments.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGARGUMENTS_HPP_
#define MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGARGUMENTS_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>

/** This is to hide away the array of arguments possible and allows
 * via getters to focus on a single entity.
 */
struct StoppingArguments
{
	/** Cstor for StoppingArguments.
	 *
	 * Sets nonsense default values
	 */
	StoppingArguments() :
		MaxWalltime(boost::chrono::seconds(0)),
		MaxIterations(0),
		Tolerance(0)
	{}

	// getter

	boost::chrono::duration<double> getMaxWalltime() const
	{ return MaxWalltime; }

	int getMaxIterations() const
	{ return MaxIterations; }

	double getTolerance() const
	{ return Tolerance; }

	// setter

	void setMaxWalltime(const boost::chrono::duration<double> _MaxWalltime)
	{ MaxWalltime = _MaxWalltime; }

	void setMaxIterations(const int _MaxIterations)
	{ MaxIterations = _MaxIterations; }

	void setTolerance(const double _Tolerance)
	{ Tolerance = _Tolerance; }

private:
	boost::chrono::duration<double> MaxWalltime;
	int MaxIterations;
	double Tolerance;
};


#endif /* MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGARGUMENTS_HPP_ */
