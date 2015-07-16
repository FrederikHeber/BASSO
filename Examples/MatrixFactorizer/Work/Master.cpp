/*
 * Master.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Master.hpp"

#ifdef MPI_FOUND
#include <boost/mpi/environment.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/nonblocking.hpp>

namespace mpi = boost::mpi;

#include <deque>

#include "Log/Logging.hpp"
#include "MatrixFactorizer/Helpers/detail.hpp"
#include "MatrixFactorizer/Work/WorkPackage.hpp"
#include "Minimizations/Elements/Eigen_matrix_serialization.hpp"
#include "Options/CommandLineOptions.hpp"

using namespace detail;

Master::Master(mpi::communicator &_world) :
	world(_world)
{}

bool Master::solve(
		const CommandLineOptions &_opts,
		const Eigen::MatrixXd &_matrix,
		const Eigen::MatrixXd &_rhs,
		Eigen::MatrixXd &_solution,
		const unsigned int _loop_nr
		)
{
	// send round that we don't terminate yet
	BOOST_LOG_TRIVIAL(debug)
			<< "#0 - broadcasting no terminate.";
	bool full_terminate = false;
	mpi::broadcast(world, full_terminate, 0);

	// gather all available workers
	std::deque<int> AvailableWorkers;
	for (int i=1;i<world.size();++i) {
		int id = -1;
		BOOST_LOG_TRIVIAL(debug)
				<< "#0 - waiting for id " << i << ".";
		world.recv(i, InitialId, id);
		assert( id != -1 );
		AvailableWorkers.push_back(id);
	}

	// send round global information
	BOOST_LOG_TRIVIAL(debug)
			<< "#0 - broadcasting matrix.";
	mpi::broadcast(world, const_cast<Eigen::MatrixXd &>(_matrix), 0);
	mpi::broadcast(world, const_cast<unsigned int &>(_loop_nr), 0);

	// allocate variables for non-blocking communication
	std::vector<WorkPackage> packages(world.size());
	std::vector<WorkResult> results(world.size());
	std::deque<mpi::request> uncompleted_requests;

	// go through all columns of rhs
	bool continue_condition = true;
	for (int col = 0;
			(continue_condition) && (col < _rhs.cols());
			++col) {
		// send current column of rhs as work package
		const size_t freeworker = AvailableWorkers.front();
		AvailableWorkers.pop_front();
		packages[freeworker] =
				WorkPackage(_rhs.col(col), _solution.col(col), col);
		BOOST_LOG_TRIVIAL(debug)
				<< "#0 - sending work package " << col
				<< " to " << freeworker << ".";
		world.isend(freeworker, ColumnWise, packages[freeworker]);

		// receive solution in a non-blocking manner
		BOOST_LOG_TRIVIAL(debug)
				<< "#0 - launching receiving work package " << col
				<< " from " << freeworker << ".";
		mpi::request request_object =
				world.irecv(freeworker, ColumnWise, results[freeworker]);
		uncompleted_requests.push_back(request_object);

		// wait for free workers if non available
		if ((AvailableWorkers.empty()) && (!uncompleted_requests.empty())) {
			// wait for any result
			BOOST_LOG_TRIVIAL(debug)
					<< "#0 - waiting for any work result.";
			std::pair<
				mpi::status,
				std::deque<mpi::request>::iterator > result =
					mpi::wait_any(
							uncompleted_requests.begin(),
							uncompleted_requests.end());
			uncompleted_requests.erase(result.second);
			BOOST_LOG_TRIVIAL(debug)
					<< "#0 - received work result from "
					<< result.first.source() << ".";
			// place worker back into available deque
			AvailableWorkers.push_back(result.first.source());
			// store solution
			continue_condition &=
					handleResult(result.first, results, _solution);
		}
	}
	// wait for all remaining results and handle incoming results
	BOOST_LOG_TRIVIAL(debug)
			<< "#0 - waiting for all remaining work results.";
	{
		std::vector<mpi::status> completed_requests;
		mpi::wait_all(
				uncompleted_requests.begin(),
				uncompleted_requests.end(),
				std::back_inserter(completed_requests));
		uncompleted_requests.clear();
		for (std::vector<mpi::status>::const_iterator iter = completed_requests.begin();
				iter != completed_requests.end(); ++iter) {
			// store solution
			continue_condition &=
					handleResult(*iter, results, _solution);
		}
	}

	// send terminate signal by empty work package
	for (int i=1;i<world.size();++i) {
		BOOST_LOG_TRIVIAL(debug)
				<< "#0 - sending termination signal to " << i << ".";
		world.send(i, ColumnWise, WorkPackage());
	}

	return continue_condition;
}

bool Master::handleResult(
		const mpi::status &_result,
		const std::vector<WorkResult> &_results,
		Eigen::MatrixXd &_solution
		)
{
	bool solve_ok = true;
	if (_result.error() != 0) {
		// look at stop_condition and store solution
		const int id = _result.source();
		solve_ok &= _results[id].solve_ok;
		_solution.col(_results[id].col) = _results[id].solution;
		if ((_solution.innerSize() > 10) || (_solution.outerSize() > 10)) {
			BOOST_LOG_TRIVIAL(trace)
					<< "Got solution for col #" << _results[id].col
					<< "\n" << _results[id].solution.transpose();
		} else {
			BOOST_LOG_TRIVIAL(debug)
					<< "Got solution for col #" << _results[id].col
					<< "\n" << _results[id].solution.transpose();
		}
	} else
		solve_ok = false;

	return solve_ok;
}

void Master::sendTerminate()
{
	// send full terminate signal
	// send round global information
	BOOST_LOG_TRIVIAL(debug)
			<< "#0 - broadcasting full terminate.";
	bool full_terminate = true;
	mpi::broadcast(world, full_terminate, 0);
}

#endif /* MPI_FOUND */

