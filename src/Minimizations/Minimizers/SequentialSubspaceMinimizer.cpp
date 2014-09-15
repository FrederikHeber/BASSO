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
#include "Minimizations/Functions/FunctionMinimizer.hpp"
#include "Minimizations/Functions/HyperplaneProjection.hpp"
#include "Minimizations/Functions/MinimizationFunctional.hpp"
#include "Minimizations/Functions/VectorProjection.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

using namespace boost::assign;

// instantiate required template functions
CONSTRUCT_FUNCTIONMINIMIZER(Eigen::VectorXd)

SequentialSubspaceMinimizer::SequentialSubspaceMinimizer(
		const InverseProblem_ptr_t &_inverseproblem,
		const double _Delta,
		const unsigned int _maxiter,
		Database &_database,
		const unsigned int _outputsteps
		) :
	GeneralMinimizer(
			_inverseproblem,
			_Delta,
			_maxiter,
			_database,
			_outputsteps
			),
	N(2),
	MatrixVectorProduct_subspace(MatrixVectorProduct),
	ScalarVectorProduct_subspace(ScalarVectorProduct),
	projector(*_inverseproblem->x->getSpace()->getDualSpace()->getNorm(),
			dynamic_cast<const PowerTypeDualityMapping &>(
					*_inverseproblem->x->getSpace()->getDualSpace()->getDualityMapping()
					),
			_inverseproblem->x->getSpace()->getDualSpace()->getDualityMapping()->getPower(),
			ScalarVectorProduct_subspace)
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
	angle_tuple.insert( std::make_pair("p", _val_NormX));
	angle_tuple.insert( std::make_pair("r", _val_NormY));
	angle_tuple.insert( std::make_pair("N", (int)_N));
	angle_tuple.insert( std::make_pair("dim", (int)_dim));
	angle_tuple.insert( std::make_pair("iteration", (int)0));
	std::vector<std::string> names;
	names += "angle","bregman_angle";
	for (std::vector<std::string>::const_iterator nameiter = names.begin();
			nameiter != names.end(); ++nameiter)
		for (unsigned int i=0; i<MAXANGLES; ++i) {
			std::stringstream componentname;
			componentname << *nameiter << i+1;
			BOOST_LOG_TRIVIAL(debug)
				<< "Adding " << componentname.str();
			angle_tuple.insert( std::make_pair(componentname.str(), 0.));
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
		istate.set(_startvalue, residual, residuum, N);
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

	// build data tuple for overall information
	Table& overall_table = database.addTable("overall");
	Table::Tuple_t overall_tuple = prepareOverallTuple(
			NormX.getPvalue(), NormY.getPvalue(), N, SpaceX.getDimension());

	// build angle tuple for search direction angle information
	Table& angle_table = database.addTable("angles");
	Table::Tuple_t angle_tuple = prepareAngleTuple(
			NormX.getPvalue(), NormY.getPvalue(), N, SpaceX.getDimension());

	/// -# check stopping criterion
	const double ynorm = NormY(y);
	bool StopCriterion = false;
	StopCriterion = (fabs(istate.residuum/ynorm) <= TolY);

	while (!StopCriterion) {
		per_iteration_tuple.replace( "iteration", (int)istate.NumberOuterIterations);
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
			const double new_error = NormX(istate.m_solution-_truesolution);
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
		// calculate bregman angles for angles database
		{
			const IterationState::angles_t angles =
					istate.calculateBregmanAngles(
							DualNormX,
							projector,
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
		// update search space with new direction
		istate.updateSearchSpace(newdir, alpha);
		per_iteration_tuple.replace( "updated_index", (int)istate.getIndex());
		BOOST_LOG_TRIVIAL(trace)
			<< "updated_index is " << istate.getIndex();

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

			HyperplaneProjection functional(
					bregman,
					dual_solution->getVectorRepresentation(),
					istate.getSearchSpace(),
					istate.getAlphas());

			FunctionMinimizer<Eigen::VectorXd> minimizer(
				functional, tmin);

			tmin = minimizer(N, TolFun, tmin);

			BOOST_LOG_TRIVIAL(trace)
				<< "tmin is " << tmin.transpose();
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

		// check iterations count
		++istate.NumberOuterIterations;
		StopCriterion =
				(istate.NumberOuterIterations >= MaxOuterIterations)
				|| (fabs(istate.residuum/ynorm) <= TolY);

		// print intermediat solution
		printIntermediateSolution(
				 istate.m_solution->getVectorRepresentation(),
				 A.getMatrixRepresentation(),
				istate.NumberOuterIterations);

		// submit current tuples
		per_iteration_table.addTuple(per_iteration_tuple);
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
		const InverseProblem_ptr_t &_problem,
		enum UpdateAlgorithmType _type)
{
	const NormedSpace & SpaceX = *_problem->A->getSourceSpace();
	const NormedSpace & DualSpaceX = *SpaceX.getDualSpace();
	const Norm & DualNormX = *DualSpaceX.getNorm();
	switch (_type) {
	case RoundRobin:
		istate.updateIndex = &IterationState::advanceIndex;
		break;
	case MostParallel:
		istate.updateIndex = boost::bind(
				&IterationState::updateIndexToMostParallel,
				_1,
				boost::cref(DualNormX),
				boost::cref(projector),
				_2);
		break;
	case MostOrthogonal:
		istate.updateIndex = boost::bind(
				&IterationState::updateIndexToMostOrthogonal,
				_1,
				boost::cref(DualNormX),
				boost::cref(projector),
				_2);
		break;
	default:
		BOOST_LOG_TRIVIAL(error)
			<< "Unknown updateIndex algorithm.";
		assert(0);
	}
}

void
SequentialSubspaceMinimizer::setEnforceRandomMapping(
		const bool _enforceRandomMapping)
{
	istate.enforceRandomMapping = _enforceRandomMapping;
}

