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

#include <cassert>
#include <deque>
#include <functional>
#include <numeric>

#include "Database/Table.hpp"
#include "Log/Logging.hpp"
#include "MatrixFactorizer/Helpers/detail.hpp"
#include "MatrixFactorizer/Work/WorkPackage.hpp"
#include "Minimizations/Elements/Eigen_matrix_serialization.hpp"
#include "Options/CommandLineOptions.hpp"

using namespace detail;

Master::Master(
		mpi::communicator &_world,
		const InnerProblemDatabase::keys_t &_overall_keys) :
	world(_world),
	overall_keys(_overall_keys),
	projector_db(new InnerProblemDatabase(_overall_keys)),
	solver_db(new InnerProblemDatabase(_overall_keys))
{
	// only master should use this class directly, others use Slave
	assert( world.rank() == 0);
}

void printNumberOfValues(
		const size_t _rank,
		const std::string &_name,
		const std::vector<AccumulatedValues> &_values
		)
{
	std::vector<size_t> sizes(_values.size(), 0);
	std::transform(
			const_cast< const std::vector<AccumulatedValues>& >(_values).begin(),
			const_cast< const std::vector<AccumulatedValues>& >(_values).end(),
			sizes.begin(),
			std::mem_fun_ref(&AccumulatedValues::getNumberOfValues));
	std::stringstream output;
	std::copy(sizes.begin(), sizes.end(),
			std::ostream_iterator<size_t>(output, ","));
	BOOST_LOG_TRIVIAL(debug)
			<< "#" << _rank << " gathered " << output.str();
	const size_t totalentries =
			std::accumulate(sizes.begin(), sizes.end(), 0);
	BOOST_LOG_TRIVIAL(debug)
			<< "#" << _rank << " gathered " << totalentries
			<< " accumulated values from " << _name << " problem.";
}

bool Master::solve(
		const CommandLineOptions &_opts,
		const Eigen::MatrixXd &_matrix,
		const Eigen::MatrixXd &_rhs,
		Eigen::MatrixXd &_solution,
		const std::string &_auxiliary_constraints
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
	assert(AvailableWorkers.size() == world.size()-(unsigned int)1);

	// send round global information
	BOOST_LOG_TRIVIAL(debug)
			<< "#0 - broadcasting options.";
	mpi::broadcast(world, const_cast<CommandLineOptions &>(_opts), 0);
	mpi::broadcast(world, const_cast<InnerProblemDatabase::keys_t &>(overall_keys), 0);
	BOOST_LOG_TRIVIAL(debug)
			<< "#0 - broadcasting constraints.";
	mpi::broadcast(world, const_cast<std::string &>(_auxiliary_constraints), 0);
	BOOST_LOG_TRIVIAL(debug)
			<< "#0 - broadcasting matrix.";
	mpi::broadcast(world, const_cast<Eigen::MatrixXd &>(_matrix), 0);

	// allocate variables for non-blocking communication
	std::vector<WorkResult> results(world.size()-1);
	std::deque<mpi::request> uncompleted_requests;

	// go through all columns of rhs
	bool solver_ok = true;
	for (int col = 0;
			(solver_ok) && (col < _rhs.cols());
			++col) {
		// send current column of rhs as work package
		const size_t freeworker = AvailableWorkers.front();
		AvailableWorkers.pop_front();
		WorkPackage package =
				WorkPackage(_rhs.col(col), _solution.col(col), col);
		BOOST_LOG_TRIVIAL(debug)
				<< "#0 - sending work package " << col
				<< " to " << freeworker << ".";
		// if we ever use isend again, then we also need to store the
		// resulting mpi::requests for isend (and not only for irecv)
		// and wait for \b both.
		// See here http://stackoverflow.com/questions/4024940/boost-mpi-whats-received-isnt-what-was-sent
		world.send(freeworker, ColumnWise, package);

		// receive solution in a non-blocking manner
		BOOST_LOG_TRIVIAL(debug)
				<< "#0 - launching receiving work package " << col
				<< " from " << freeworker << ".";
		mpi::request request_object =
				world.irecv(freeworker, ColumnWise, results[freeworker-1]);
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
			solver_ok &=
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
			solver_ok &=
					handleResult(*iter, results, _solution);
		}
	}

	// send terminate signal by empty work package
	for (int i=1;i<world.size();++i) {
		BOOST_LOG_TRIVIAL(debug)
				<< "#0 - sending termination signal to " << i << ".";
		world.send(i, ColumnWise, WorkPackage());
	}
	// gather operation for values accumulated from inner problem's overall tables
	{
		std::vector<AccumulatedValues> out_values;
		out_values.resize(world.size());
		AccumulatedValues in_values;
		boost::mpi::gather(world, in_values, out_values, 0);
		printNumberOfValues(world.rank(), "projector", out_values);
		// skip the master's empty values
		static_cast<InnerProblemDatabase &>(*projector_db.get()).
				insertValues(++out_values.begin(), out_values.end());
	}
	{
		std::vector<AccumulatedValues> out_values;
		out_values.resize(world.size());
		AccumulatedValues in_values;
		boost::mpi::gather(world, in_values, out_values, 0);
		printNumberOfValues(world.rank(), "minimization", out_values);
		// skip the master's empty values
		static_cast<InnerProblemDatabase &>(*solver_db.get()).
				insertValues(++out_values.begin(), out_values.end());
	}

	return solver_ok;
}

void Master::insertAccumulatedProjectorValues(
		Table &_table,
		const std::string &_suffix) const
{
	static_cast<InnerProblemDatabase &>(*projector_db.get()).
			insertAccumulatedValues(_table, _suffix);
}

void Master::insertAccumulatedSolverValues(
		Table &_table,
		const std::string &_suffix) const
{
	static_cast<InnerProblemDatabase &>(*solver_db.get()).
			insertAccumulatedValues(_table, _suffix);
}

void Master::resetAccumulatedProjectorValues()
{
	projector_db.reset(new InnerProblemDatabase(overall_keys));
}

void Master::resetAccumulatedSolverValues()
{
	solver_db.reset(new InnerProblemDatabase(overall_keys));
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
		solve_ok &= _results[id-1].solve_ok;
		_solution.col(_results[id-1].col) = _results[id-1].solution;
		if ((_solution.innerSize() > 10) || (_solution.outerSize() > 10)) {
			BOOST_LOG_TRIVIAL(trace)
					<< "Got solution for col #" << _results[id-1].col
					<< "\n" << _results[id-1].solution.transpose();
		} else {
			BOOST_LOG_TRIVIAL(debug)
					<< "Got solution for col #" << _results[id-1].col
					<< "\n" << _results[id-1].solution.transpose();
		}
	} else {
		BOOST_LOG_TRIVIAL(error)
				<< "Could not open solution for column "
				<< _results[_result.source()-1].col;
		solve_ok = false;
	}

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

