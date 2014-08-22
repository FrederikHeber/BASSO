/*
 * OperationCounter.hpp
 *
 *  Created on: Aug 22, 2014
 *      Author: heber
 */

#ifndef OPERATIONCOUNTER_HPP_
#define OPERATIONCOUNTER_HPP_

#include <boost/bind.hpp>
#include <boost/chrono.hpp>
#include <boost/function.hpp>

template <class res, class T, class V>
class OperationCounter
{
public:
	//!> typedef for the bound function
	typedef boost::function<res(T, V)> function_t;

	OperationCounter(function_t &_operatorfunction) :
		operatorfunction(_operatorfunction),
		count(0),
		timing(boost::chrono::nanoseconds(0))
	{}
	~OperationCounter()
	{}

	res operator()(T _a, V _b) const
	{
		++count;
		boost::chrono::high_resolution_clock::time_point timing_start =
				boost::chrono::high_resolution_clock::now();
		res tmp = operatorfunction(_a,_b);
		boost::chrono::high_resolution_clock::time_point timing_end =
				boost::chrono::high_resolution_clock::now();
		timing += timing_end - timing_start;
		return tmp;
	}

private:
	//!> function call that is wrapped by this class
	function_t operatorfunction;

	//!> counts the number of times the operation has been called.
	mutable size_t count;

	//!> counts the time used for this operation
	mutable boost::chrono::nanoseconds timing;
};

#endif /* OPERATIONCOUNTER_HPP_ */
