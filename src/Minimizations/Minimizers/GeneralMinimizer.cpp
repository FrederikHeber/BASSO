/*
 * GeneralMinimizer.cpp
 *
 *  Created on: Apr 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "GeneralMinimizer.hpp"

#include "MatrixIO/MatrixIO.hpp"
#include "MatrixIO/OperationCounter.hpp"

#include <boost/assign.hpp>
#include <boost/log/trivial.hpp>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/NormFactory.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMappingFactory.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

using namespace boost::assign;

// static entities
GeneralMinimizer::MinLib_names_t GeneralMinimizer::MinLib_names;

GeneralMinimizer::GeneralMinimizer(
		const InverseProblem_ptr_t &_inverseproblem,
		const double _Delta,
		const unsigned int _maxiter,
		const unsigned int _maxinneriter,
		Database &_database,
		const unsigned int _outputsteps
		) :
	Delta(_Delta),
	MaxWalltime(0.),
	MaxOuterIterations(_maxiter),
	MaxInnerIterations(_maxinneriter),
	TolX(1e-6),
	TolY(Delta),
	TolFun(1e-12),
	outputsteps(_outputsteps),
	MinLib(gnuscientificlibrary),
	database(_database)
{
	// set tolerances values
	_inverseproblem->x->getSpace()->getDualityMapping()->setTolerance(TolX);
	_inverseproblem->x->getSpace()->getDualSpace()->getDualityMapping()->setTolerance(TolX);
	_inverseproblem->y->getSpace()->getDualityMapping()->setTolerance(TolY);

	// initalize list with static names
	if (MinLib_names.empty())
		MinLib_names +=
			std::make_pair( "gsl", gnuscientificlibrary),
			std::make_pair( "nlopt", nonlinearoptimization);

}

bool GeneralMinimizer::CheckWalltime(
		const boost::chrono::duration<double> &_time) const
{
	if (MaxWalltime != boost::chrono::duration<double>(0.))
		return _time >= MaxWalltime;
	else
		return false;
}

bool GeneralMinimizer::CheckIterations(
		const int _current_outeriterations) const
{
	// walltime overrules maxiter
	if (MaxWalltime == boost::chrono::duration<double>(0.))
		return _current_outeriterations >= MaxOuterIterations;
	else
		return false;
}

bool GeneralMinimizer::CheckResiduum(
		const double _residuum) const
{
	return _residuum <= TolY;
}

double GeneralMinimizer::calculateResidual(
		const InverseProblem_ptr_t &_problem,
		SpaceElement_ptr_t &_residual
		) const
{
	const LinearMapping &A = static_cast<const LinearMapping &>(*_problem->A);
	*_residual = A * _problem->x;
	*_residual -= _problem->y;
	const Norm &NormY = *_problem->y->getSpace()->getNorm();
	return NormY(_residual);
}

bool GeneralMinimizer::isValidMinLibName(const std::string &_name)
{
	MinLib_names_t::const_iterator iter =
			MinLib_names.find(_name);
	return (iter != MinLib_names.end());
}

void GeneralMinimizer::setMinLib(const std::string &_name)
{
	assert( isValidMinLibName(_name) );
	MinLib = MinLib_names[_name];
}

void GeneralMinimizer::printIntermediateSolution(
		const SpaceElement_ptr_t &_solution,
		const LinearMapping &_A,
		unsigned int _NumberOuterIterations
		) const
{
	// print each solution
	if ((outputsteps != 0) &&
			(_NumberOuterIterations % outputsteps == 0)) {
		{
			std::stringstream solution_file;
			solution_file << "solution"
					<< (_NumberOuterIterations / outputsteps) << ".m";
			using namespace MatrixIO;
			std::ofstream ost(solution_file.str().c_str());
			if (ost.good())
				try {
					ost << _solution;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Could not write all data of intermediate solution to stream.\n";
				}
			else {
				std::cerr << "Failed to open " << solution_file.str() << std::endl;
			}
		}
		{
			std::stringstream solution_file;
			solution_file << "projected_solution"
					<< (_NumberOuterIterations / outputsteps) << ".m";
			using namespace MatrixIO;
			std::ofstream ost(solution_file.str().c_str());
			if (ost.good())
				try {
					ost << _A * _solution;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Could not write all data of projected intermediate solution to stream.\n";
				}
			else {
				std::cerr << "Failed to open " << solution_file.str() << std::endl;
			}
		}

	}
}

Table::Tuple_t GeneralMinimizer::preparePerIterationTuple(
		const double _val_NormX,
		const double _val_NormY,
		const unsigned int _N,
		const unsigned int _dim,
		const int _MaxOuterIterations)
{
	Table::Tuple_t per_iteration_tuple;
	per_iteration_tuple.insert( std::make_pair("p", _val_NormX), Table::Parameter);
	per_iteration_tuple.insert( std::make_pair("r", _val_NormY), Table::Parameter);
	per_iteration_tuple.insert( std::make_pair("N", (int)_N), Table::Parameter);
	per_iteration_tuple.insert( std::make_pair("dim", (int)_dim), Table::Parameter);
	per_iteration_tuple.insert( std::make_pair("max_iterations", _MaxOuterIterations), Table::Parameter);
	per_iteration_tuple.insert( std::make_pair("iteration", (int)0), Table::Data);
	per_iteration_tuple.insert( std::make_pair("stepwidth", 0.), Table::Data);
	per_iteration_tuple.insert( std::make_pair("relative_residual", 0.), Table::Data);
	per_iteration_tuple.insert( std::make_pair("error", 0.), Table::Data);
	per_iteration_tuple.insert( std::make_pair("bregman_distance", 0.), Table::Data);
	per_iteration_tuple.insert( std::make_pair("updated_index", (int)0), Table::Data);
	return per_iteration_tuple;
}

Table::Tuple_t GeneralMinimizer::prepareOverallTuple(
		const double _val_NormX,
		const double _val_NormY,
		const unsigned int _N,
		const unsigned int _dim,
		const int _MaxOuterIterations)
{
	Table::Tuple_t overall_tuple;
	overall_tuple.insert( std::make_pair("p", _val_NormX), Table::Parameter);
	overall_tuple.insert( std::make_pair("r", _val_NormY), Table::Parameter);
	overall_tuple.insert( std::make_pair("N", (int)_N), Table::Parameter);
	overall_tuple.insert( std::make_pair("dim", (int)_dim), Table::Parameter);
	overall_tuple.insert( std::make_pair("max_iterations", _MaxOuterIterations), Table::Parameter);
	overall_tuple.insert( std::make_pair("iterations", (int)0), Table::Data);
	overall_tuple.insert( std::make_pair("relative_residual", 0.), Table::Data);
	overall_tuple.insert( std::make_pair("runtime", 0.), Table::Data);
	overall_tuple.insert( std::make_pair("element_creation_operations", (int)0), Table::Data);
	overall_tuple.insert( std::make_pair("linear_time_operations", (int)0), Table::Data);
	overall_tuple.insert( std::make_pair("quadratic_time_operations", (int)0), Table::Data);
	overall_tuple.insert( std::make_pair("element_creation_runtime", 0.), Table::Data);
	overall_tuple.insert( std::make_pair("linear_time_runtime", 0.), Table::Data);
	overall_tuple.insert( std::make_pair("quadratic_time_runtime", 0.), Table::Data);
	return overall_tuple;
}

