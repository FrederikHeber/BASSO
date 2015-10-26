/*
 * SolverFactory.cpp
 *
 *  Created on: Jan 28, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "SolverFactory.hpp"

#include <boost/assign.hpp>

#include "Options/CommandLineOptions.hpp"
#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/InverseProblems/InverseProblemFactory.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/LandweberMinimizer.hpp"
#include "Minimizations/Minimizers/SequentialSubspaceMinimizer.hpp"
#include "Minimizations/Minimizers/SequentialSubspaceMinimizerNoise.hpp"
#include "Minimizations/Minimizers/StepWidths/DetermineStepWidthFactory.hpp"

using namespace boost::assign;

InverseProblem_ptr_t SolverFactory::createInverseProblem(
		const CommandLineOptions &_opts,
		const Eigen::MatrixXd &_matrix,
		const Eigen::VectorXd &_rhs)
{
	InverseProblem_ptr_t inverseproblem;
	InverseProblemFactory::args_t args_SpaceX;
	InverseProblemFactory::args_t args_SpaceY;
	// starts with "lp.."
	if (_opts.type_spacex.find("lp") == 0) {
		args_SpaceX +=
				boost::any(_opts.px),
				boost::any(_opts.powerx);
	} else if (_opts.type_spacex == "regularized_l1") {
		args_SpaceX +=
				boost::any(_opts.regularization_parameter),
				boost::any(_opts.powerx);
	}
	// starts with "lp.."
	if (_opts.type_spacey.find("lp") == 0) {
		args_SpaceY +=
				boost::any(_opts.py),
				boost::any(_opts.powery);
	}
	inverseproblem = InverseProblemFactory::create(
			_opts.type_spacex,
			args_SpaceX,
			_opts.type_spacey,
			args_SpaceY,
			_matrix,
			_rhs);

	return inverseproblem;
}

Database_ptr_t SolverFactory::createDatabase(
		const CommandLineOptions &_opts
		)
{
	Database_ptr_t database(new Database());
	if (!_opts.iteration_file.string().empty())
		database->setDatabaseFile(_opts.iteration_file.string());
	database->setReplacePresentParameterTuples(_opts.database_replace);

	return database;
}

MinimizerFactory::instance_ptr_t SolverFactory::createMinimizer(
		const CommandLineOptions &_opts,
		InverseProblem_ptr_t &_inverseproblem,
		Database_ptr_t &_database,
		const StoppingCriterion::ptr_t &_stopping_criteria,
		const unsigned int _maxiterations
		)
{
	MinimizerFactory factory;
	MinimizerFactory::instance_ptr_t minimizer;
	// set regularization parameter in case of regularizedl1norm
	minimizer =
		factory.createInstance(
			_opts.type,
			_inverseproblem,
			_opts.delta,
			_maxiterations,
			_opts.maxinneriter,
			*_database,
			_stopping_criteria,
			(const enum DetermineStepWidthFactory::stepwidth_enumeration)_opts.stepwidth_type,
			_opts.outputsteps,
			_opts.orthogonalization_type);
	minimizer->setMinLib(_opts.minlib);
	// hand over additional parameters to add to submitted tuples
	minimizer->setAdditionalTupleParameters(_opts.tuple_parameters);

	try {
		// create instance with some specifics
		switch(_opts.type) {
		case MinimizerFactory::landweber:
			static_cast<LandweberMinimizer*>(
					minimizer.get())->setC(_opts.C);
			break;
		case MinimizerFactory::sequentialsubspace:
			static_cast<SequentialSubspaceMinimizer*>(minimizer.get())->setN(_opts.N);
			static_cast<SequentialSubspaceMinimizer*>(
					minimizer.get())->setupdateIndexAlgorithm(
							_opts.updatetype);
			static_cast<SequentialSubspaceMinimizer*>(
					minimizer.get())->setEnforceRandomMapping(
							_opts.enforceRandomMapping);
			static_cast<SequentialSubspaceMinimizer*>(minimizer.get())->setInexactLinesearch(
					_opts.inexactLinesearch);
			if (_opts.wolfe_constants.size() == 2) {
				static_cast<SequentialSubspaceMinimizer*>(minimizer.get())->setWolfeConstants(
						_opts.wolfe_constants);
				// warning in case sanity check fails
				if (!_opts.inexactLinesearch)
					BOOST_LOG_TRIVIAL(warning)
						<< "Wolfe constants set although we do perform an exact line search.";
			}
			static_cast<SequentialSubspaceMinimizer*>(
					minimizer.get())->setDoCalculateAngles(
							_opts.calculateAngles);
			break;
		case MinimizerFactory::sequentialsubspace_noise:
			static_cast<SequentialSubspaceMinimizerNoise*>(
					minimizer.get())->setN(_opts.N);
			break;
		default:
			std::cerr << "Unknown InstanceType"
				<< MinimizerFactory::getNameForType(_opts.type) << "." << std::endl;
			return MinimizerFactory::instance_ptr_t();
			break;
		}
	} catch (MinimizationIllegalValue_exception &e) {
		std::cerr << "Illegal value for "
				<< *boost::get_error_info<MinimizationIllegalValue_name>(e)
				<< std::endl;
		return MinimizerFactory::instance_ptr_t();
	}

	return minimizer;
}
