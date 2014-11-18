/*
 * SequentialSubspaceMinimizer.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "SequentialSubspaceMinimizer.hpp"

#include <boost/chrono.hpp>
#include <boost/assign.hpp>
#include <boost/log/trivial.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>

#include "Database/Database.hpp"
#include "Database/Table.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/Functions/HyperplaneProjection.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizer.hpp"
#include "Minimizations/Functions/Minimizers/MinimizationFunctional.hpp"
#include "Minimizations/Functions/Minimizers/Minimizer.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/NormFactory.hpp"
#include "Minimizations/Norms/RegularizedL1Norm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

using namespace boost::assign;

SequentialSubspaceMinimizer::SequentialSubspaceMinimizer(
		const InverseProblem_ptr_t &_inverseproblem,
		const double _Delta,
		const unsigned int _maxiter,
		const unsigned int _maxinneriter,
		Database &_database,
		const unsigned int _outputsteps
		) :
	GeneralMinimizer(
			_inverseproblem,
			_Delta,
			_maxiter,
			_maxinneriter,
			_database,
			_outputsteps
			),
	N(2),
	MatrixVectorProduct_subspace(MatrixVectorProduct),
	ScalarVectorProduct_subspace(ScalarVectorProduct),
	inexactLinesearch(false),
	constant_positivity(1e-6),
	constant_interpolation(0.6),
	DoCalculateAngles(false)
{}

void SequentialSubspaceMinimizer::setN(
		const unsigned int _N
		)
{
	// check that number of search directions is greater than 1
	if (_N < 1)
		throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("N");
	const_cast<unsigned int&>(N) = _N;
	// make state invalid
	istate.reset();
}

//!> limit for angle columns in angles table
static const unsigned int MAXANGLES = 4;

static Table::Tuple_t prepareAngleTuple(
		const double _val_NormX,
		const double _val_NormY,
		const unsigned int _N,
		const unsigned int _dim)
{
	Table::Tuple_t angle_tuple;
	angle_tuple.insert( std::make_pair("p", _val_NormX), Table::Parameter);
	angle_tuple.insert( std::make_pair("r", _val_NormY), Table::Parameter);
	angle_tuple.insert( std::make_pair("N", (int)_N), Table::Parameter);
	angle_tuple.insert( std::make_pair("dim", (int)_dim), Table::Parameter);
	angle_tuple.insert( std::make_pair("iteration", (int)0), Table::Data);
	std::vector<std::string> names;
	names += "angle","bregman_angle";
	for (std::vector<std::string>::const_iterator nameiter = names.begin();
			nameiter != names.end(); ++nameiter)
		for (unsigned int i=0; i<MAXANGLES; ++i) {
			std::stringstream componentname;
			componentname << *nameiter << i+1;
			BOOST_LOG_TRIVIAL(debug)
				<< "Adding " << componentname.str();
			angle_tuple.insert( std::make_pair(componentname.str(), 0.), Table::Data);
		}
	return angle_tuple;
}

SequentialSubspaceMinimizer::ReturnValues
SequentialSubspaceMinimizer::operator()(
		const InverseProblem_ptr_t &_problem,
		const SpaceElement_ptr_t &_startvalue,
		const SpaceElement_ptr_t &_dualstartvalue,
		const SpaceElement_ptr_t &_truesolution
		)
{
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	// gather some refs for convenient access
	const NormedSpace & SpaceX = *_problem->A->getSourceSpace();
	const NormedSpace & DualSpaceX = *SpaceX.getDualSpace();
	const NormedSpace & SpaceY = *_problem->A->getTargetSpace();
	const Norm & NormX = *SpaceX.getNorm();
	const Norm & DualNormX = *DualSpaceX.getNorm();
	const Norm & NormY = *SpaceY.getNorm();
	const Mapping & J_p = *SpaceX.getDualityMapping();
	const Mapping & J_q = *DualSpaceX.getDualityMapping();
	const Mapping & j_r = *SpaceY.getDualityMapping();
	const SpaceElement_ptr_t &y = _problem->y;
	const LinearMapping &A =
			dynamic_cast<const LinearMapping &>(*_problem->A);
	const Mapping_ptr_t A_adjoint = A.getAdjointMapping();
	const LinearMapping &A_t =
			dynamic_cast<const LinearMapping &>(*A_adjoint);

	/// -# initialize return structure
	if (!istate.getisInitialized()) {
		SpaceElement_ptr_t residual = A.getTargetSpace()->createElement();
		const double residuum = calculateResidual(
			_problem,
			residual);
		istate.set(
				_startvalue,
				residual,
				residuum,
				N,
				MatrixVectorProduct_subspace,
				ScalarVectorProduct_subspace);
	}

	/// -# calculate some values prior to loop
	// Jx=DualityMapping(x,NormX,PowerX,TolX);
	SpaceElement_ptr_t dual_solution = DualSpaceX.createElement();
	*dual_solution = _dualstartvalue;
	BOOST_LOG_TRIVIAL(trace)
		<< "Jx_0 is " << dual_solution;

	if ((_truesolution->getSpace()->getDimension() > 10)) {
		BOOST_LOG_TRIVIAL(trace)
			<< "true solution is " << _truesolution << std::endl;
	} else {
		BOOST_LOG_TRIVIAL(info)
			<< "true solution is " << _truesolution << std::endl;
	}

	BregmanDistance Delta_p(
			NormX,
			dynamic_cast<const PowerTypeDualityMapping &>(J_p),
			J_p.getPower(),
			ScalarVectorProduct);
	double old_distance = 0.;
	if (!_truesolution->isZero()) {
		old_distance = Delta_p(
			istate.m_solution,
			_truesolution,
			dual_solution)
			+ 1e4*BASSOTOLERANCE; // make sure its larger
		BOOST_LOG_TRIVIAL(debug)
				<< "Starting Bregman distance is " << old_distance;
	}

	// build data tuple for iteration information
	Table& per_iteration_table = database.addTable("per_iteration");
	Table::Tuple_t per_iteration_tuple = preparePerIterationTuple(
			NormX.getPvalue(), NormY.getPvalue(), N, SpaceX.getDimension());
	if (inexactLinesearch) {
		per_iteration_tuple.insert( std::make_pair("c1", constant_positivity), Table::Parameter);
		per_iteration_tuple.insert( std::make_pair("c2", constant_interpolation), Table::Parameter);
	}
	per_iteration_tuple.insert( std::make_pair("inner_iterations", (int)0), Table::Data);

	// build data tuple for overall information
	Table& overall_table = database.addTable("overall");
	Table::Tuple_t overall_tuple = prepareOverallTuple(
			NormX.getPvalue(), NormY.getPvalue(), N, SpaceX.getDimension());
	if (inexactLinesearch) {
		overall_tuple.insert( std::make_pair("c1", constant_positivity), Table::Parameter);
		overall_tuple.insert( std::make_pair("c2", constant_interpolation), Table::Parameter);
	}
	overall_tuple.insert( std::make_pair("matrix_vector_products_subspace", (int)0), Table::Data );
	overall_tuple.insert( std::make_pair("vector_vector_products_subspace", (int)0), Table::Data );

	Table& angle_table = database.addTable("angles");
	Table::Tuple_t angle_tuple;
	if (DoCalculateAngles) {
		// build angle tuple for search direction angle information
		angle_tuple = prepareAngleTuple(
				NormX.getPvalue(), NormY.getPvalue(), N, SpaceX.getDimension());
	}

	/// -# check stopping criterion
	const double ynorm = NormY(y);
	bool StopCriterion = false;
	const double initial_relative_residuum = fabs(istate.residuum/ynorm);
	StopCriterion = (initial_relative_residuum <= TolY);

	// create L2 norm for measuring error
	const Norm_ptr_t l2norm = NormFactory::createLpInstance(
			_problem->A->getSourceSpace(), 2.);

	while (!StopCriterion) {
		per_iteration_tuple.replace( "iteration", (int)istate.NumberOuterIterations);
		if (DoCalculateAngles)
			angle_tuple.replace( "iteration", (int)istate.NumberOuterIterations);
		BOOST_LOG_TRIVIAL(debug)
				<< "#" << istate.NumberOuterIterations
				<< " with residual of " << istate.residuum;
		BOOST_LOG_TRIVIAL(debug)
			<< "#" << istate.NumberOuterIterations << ": "
			<< "||Ax_n-y||/||y|| is " << istate.residuum/ynorm;
		per_iteration_tuple.replace( "relative_residual", istate.residuum/ynorm);
		// check that distance truly decreases
		if (!_truesolution->isZero()) {
			const double new_distance =
					Delta_p(
							istate.m_solution,
							_truesolution,
							dual_solution);
			BOOST_LOG_TRIVIAL(debug)
				<< "#" << istate.NumberOuterIterations << ": "
				<< "Delta_p^{x^*_n}(x_n,x) is "
				<< new_distance;
			per_iteration_tuple.replace( "bregman_distance", new_distance);
//			assert( old_distance > new_distance );
			old_distance = new_distance;
			const double new_error =
					static_cast<const RegularizedL1Norm *>(&NormX) == NULL ?
					NormX(istate.m_solution-_truesolution) :
					(*l2norm)(istate.m_solution-_truesolution);
			BOOST_LOG_TRIVIAL(debug)
				<< "#" << istate.NumberOuterIterations << ": "
				<< "||x_n-x|| is " << new_error;
			per_iteration_tuple.replace( "error", new_error);
		}
		BOOST_LOG_TRIVIAL(trace)
				<< "x_n is " << istate.m_solution;
		BOOST_LOG_TRIVIAL(trace)
				<< "R_n is " << istate.m_residual;

		// Jw=DualityMapping(w,NormY,PowerY,TolX);
		const SpaceElement_ptr_t Jw = j_r( istate.m_residual );
		BOOST_LOG_TRIVIAL(trace)
			<< "Jw= j_r (R_n) is " << Jw;

		// JwNorm=norm(w,DualNormX);
//		const double JwNorm = DualNormX(Jw);
//		BOOST_LOG_TRIVIAL(trace)
//			<< "wNorm is " << wNorm;

		// alpha=Jw'*y
		const double alpha =
				Jw * y;
		BOOST_LOG_TRIVIAL(trace)
			<< "alpha is " << alpha;

		// add u to U and alpha to alphas
		SpaceElement_ptr_t newdir =
				dual_solution->getSpace()->createElement();
		*newdir =
				MatrixVectorProduct(
						A_t.getMatrixRepresentation(),
						Jw->getVectorRepresentation());
		if (DoCalculateAngles) {
			// calculate bregman angles for angles database
			{
				const IterationState::angles_t angles =
						istate.calculateBregmanAngles(
								DualNormX,
								newdir);
				for (unsigned int i=0; (i<MAXANGLES) && (i<angles.size()); ++i) {
					std::stringstream componentname;
					componentname << "bregman_angle" << i+1;
					angle_tuple.replace( componentname.str(), angles[i]);
				}
			}
			// calculate "scalar product" angles for angles database
			{
				const IterationState::angles_t angles =
						istate.calculateAngles(
								DualNormX,
								newdir);
				for (unsigned int i=0; (i<MAXANGLES) && (i<angles.size()); ++i) {
					std::stringstream componentname;
					componentname << "angle" << i+1;
					angle_tuple.replace( componentname.str(), angles[i]);
				}
			}
		}
		// update search space with new direction
		if (_truesolution->isZero()) {
			istate.updateSearchSpace(
					newdir,
					alpha,
					dual_solution,
					istate.m_solution);
		} else {
			istate.updateSearchSpace(
					newdir,
					alpha,
					dual_solution,
					_truesolution);
		}
		per_iteration_tuple.replace( "updated_index", (int)istate.searchspace->getIndex());
		BOOST_LOG_TRIVIAL(trace)
			<< "updated_index is " << istate.searchspace->getIndex();

		Eigen::VectorXd tmin(N);
		tmin.setZero();
		{
			// tmin=fminunc(@(t) BregmanProjectionFunctional(t,Jx,u,alpha+d,DualNormX,DualPowerX,TolX),t0,BregmanOptions);
			BregmanProjectionFunctional bregman(
					DualNormX,
					dynamic_cast<const PowerTypeDualityMapping &>(J_q),
					J_q.getPower(),
					MatrixVectorProduct_subspace,
					ScalarVectorProduct_subspace);

			const HyperplaneProjection functional(
					bregman,
					dual_solution->getVectorRepresentation(),
					istate.getSearchSpace(),
					istate.getAlphas());

			// due to templation we need to instantiate both, as user
			// decides during runtime which we need
			Minimizer<gsl_vector> minimizer_gsl(N);
			Minimizer<NLopt_vector> minimizer_nlopt(N);
			const unsigned int current_index =
					istate.searchspace->getIndex();
			Wolfe_indexset_t Wolfe_indexset(1, current_index);
			if (MinLib == gnuscientificlibrary) {
				minimizer_gsl.setMaxIterations(MaxInnerIterations);
			} else if (MinLib == nonlinearoptimization) {
				// hand index set to minimizer
				minimizer_nlopt.setConstantPositivity(constant_positivity);
				minimizer_nlopt.setPositivityBoundIndices(Wolfe_indexset);

				minimizer_nlopt.setMaxIterations(MaxInnerIterations);
			} 

			unsigned int inner_iterations = 0;
			if (inexactLinesearch) {
				FunctionalMinimizer<Eigen::VectorXd, gsl_vector> fmin_gsl(
					functional,
					minimizer_gsl,
					tmin,
					constant_positivity,
					constant_interpolation);
				FunctionalMinimizer<Eigen::VectorXd, NLopt_vector> fmin_nlopt(
					functional,
					minimizer_nlopt,
					tmin,
					constant_positivity,
					constant_interpolation);

				switch (MinLib) {
				case gnuscientificlibrary:
					inner_iterations = fmin_gsl(
							N,
							TolFun,
							Wolfe_indexset,
							tmin
							);
					break;
				case nonlinearoptimization:
					inner_iterations = fmin_nlopt(
							N,
							TolFun,
							Wolfe_indexset,
							tmin
							);
					break;
				default:
					throw MinimizationIllegalValue_exception()
							<< MinimizationIllegalValue_name("MinLib");
					break;
				}

			} else {
				FunctionalMinimizer<Eigen::VectorXd, gsl_vector> fmin_gsl(
						functional, minimizer_gsl, tmin);
				FunctionalMinimizer<Eigen::VectorXd, NLopt_vector> fmin_nlopt(
						functional, minimizer_nlopt, tmin);

				switch (MinLib) {
				case gnuscientificlibrary:
					inner_iterations = fmin_gsl(
							N,
							TolFun,
							tmin);
					break;
				case nonlinearoptimization:
					inner_iterations = fmin_nlopt(
							N,
							TolFun,
							tmin);
					break;
				default:
					throw MinimizationIllegalValue_exception()
							<< MinimizationIllegalValue_name("MinLib");
					break;
				}
			}
			per_iteration_tuple.replace( "inner_iterations", (int)inner_iterations);

			BOOST_LOG_TRIVIAL(debug)
			<< "tmin " << tmin.transpose()
			<< " found in " << inner_iterations << " iterations.";
		}
		per_iteration_tuple.replace( "stepwidth", tmin.norm());
		// x=DualityMapping(Jx-tmin*u,DualNormX,DualPowerX,TolX);
		{
			const SpaceElement_ptr_t tempelement = DualSpaceX.createElement();
			*tempelement = MatrixVectorProduct_subspace(istate.getSearchSpace(),tmin);
			*dual_solution -= tempelement;
		}
		BOOST_LOG_TRIVIAL(trace)
				<< "x^*_n+1 is " << dual_solution;
		*istate.m_solution = J_q(dual_solution);
		BOOST_LOG_TRIVIAL(trace)
				<< "x_n+1 is " << istate.m_solution;
		*_problem->x = istate.m_solution;

		// update residual
		istate.residuum = calculateResidual(
					_problem,
						istate.m_residual);

		// check iterations count/wall time
		boost::chrono::high_resolution_clock::time_point timing_intermediate =
				boost::chrono::high_resolution_clock::now();
		++istate.NumberOuterIterations;
		const double current_relative_residuum = fabs(istate.residuum/ynorm);
		StopCriterion =
				CheckIterations(istate.NumberOuterIterations)
				|| CheckResiduum(current_relative_residuum)
				|| CheckWalltime(boost::chrono::duration<double>(timing_intermediate - timing_start));

		// check for non-convergence
		if (current_relative_residuum >= 1e4) {
			BOOST_LOG_TRIVIAL(info)
					<< "Current relative residuum is "
					<< current_relative_residuum
					<< " exceeding initial value of 1. by "
					<< current_relative_residuum/initial_relative_residuum;
			BOOST_LOG_TRIVIAL(error)
					<< "STOPPING ITERATION";

			// create a sequence of entries in Database till MaxOuterIterations
			for (;istate.NumberOuterIterations < MaxOuterIterations;
					++istate.NumberOuterIterations) {
				// update iterations in tuple
				per_iteration_tuple.replace( "iteration", (int)istate.NumberOuterIterations);
				if (DoCalculateAngles)
					angle_tuple.replace( "iteration", (int)istate.NumberOuterIterations);

				// submit current tuples
				per_iteration_table.addTuple(per_iteration_tuple);
				if (DoCalculateAngles)
					angle_table.addTuple(angle_tuple);
			}
		}

		// print intermediate solution
		printIntermediateSolution(
				 istate.m_solution->getVectorRepresentation(),
				 A.getMatrixRepresentation(),
				istate.NumberOuterIterations);

		// submit current tuples
		per_iteration_table.addTuple(per_iteration_tuple);
		if (DoCalculateAngles)
			angle_table.addTuple(angle_tuple);
	}

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	std::cout << "The operation took " << boost::chrono::duration<double>(timing_end - timing_start) << "." << std::endl;

	// submit overall_tuple
	overall_tuple.replace( "iterations", istate.NumberOuterIterations );
	overall_tuple.replace( "relative_residual", istate.residuum );
	overall_tuple.replace( "runtime",
			boost::chrono::duration_cast<boost::chrono::duration<double> >(timing_end - timing_start).count() );
	overall_tuple.replace( "matrix_vector_products", (int)MatrixVectorProduct.getCount() );
	overall_tuple.replace( "vector_vector_products", (int)ScalarVectorProduct.getCount() );
	overall_tuple.replace( "matrix_vector_products_subspace", (int)MatrixVectorProduct_subspace.getCount() );
	overall_tuple.replace( "vector_vector_products_subspace", (int)ScalarVectorProduct_subspace.getCount() );
	overall_table.addTuple(overall_tuple);

	// and return solution
	return static_cast<ReturnValues &>(istate);
}

void
SequentialSubspaceMinimizer::setupdateIndexAlgorithm(
		enum LastNSearchDirections::UpdateAlgorithmType _type)
{
	LastNSearchDirections::updateIndexType = _type;
}

void
SequentialSubspaceMinimizer::setEnforceRandomMapping(
		const bool _enforceRandomMapping)
{
	LastNSearchDirections::enforceRandomMapping = _enforceRandomMapping;
}

