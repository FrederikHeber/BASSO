/*
 * SolutionFactory.cpp
 *
 *  Created on: Jan 28, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "SolutionFactory.hpp"

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

InverseProblem_ptr_t SolutionFactory::createInverseProblem(
		const CommandLineOptions &_opts,
		const Eigen::MatrixXd &_matrix,
		const Eigen::VectorXd &_rhs)
{
	// prepare inverse problem
	InverseProblem_ptr_t inverseproblem;
	switch (_opts.dualitytype) {
	case CommandLineOptions::regularizedl1norm:
		inverseproblem = InverseProblemFactory::createRegularizedL1Instance(
				_opts.regularization_parameter,
				_opts.powerx,
				_opts.normy,
				_opts.powery,
				_matrix,
				_rhs);
		break;
	case CommandLineOptions::defaulttype:
	default:
		inverseproblem = InverseProblemFactory::createLpInstance(
				_opts.normx,
				_opts.powerx,
				_opts.normy,
				_opts.powery,
				_matrix,
				_rhs);
		break;
	}

	return inverseproblem;
}

Database_ptr_t SolutionFactory::createDatabase(
		const CommandLineOptions &_opts
		)
{
	Database_ptr_t database(new Database());
	if (!_opts.iteration_file.string().empty())
		database->setDatabaseFile(_opts.iteration_file.string());
	database->setReplacePresentParameterTuples(_opts.database_replace);

	return database;
}

MinimizerFactory::instance_ptr_t SolutionFactory::createMinimizer(
		const CommandLineOptions &_opts,
		InverseProblem_ptr_t &_inverseproblem,
		Database_ptr_t &_database,
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
					minimizer.get())->setTau(_opts.tau);
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
