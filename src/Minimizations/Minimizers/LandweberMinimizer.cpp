/*
 * LandweberMinimizer.cpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "LandweberMinimizer.hpp"

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/chrono.hpp>
#include <boost/log/trivial.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <limits>
#include <fstream>
#include <sstream>

#include "Database/Database.hpp"
#include "Database/Table.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/Functions/ResidualFunctional.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Minimizers/StepWidths/DetermineStepWidth.hpp"
#include "Minimizations/Minimizers/StepWidths/DetermineStepWidthFactory.hpp"
#include "Minimizations/Norms/Norm.hpp"

LandweberMinimizer::LandweberMinimizer(
		const InverseProblem_ptr_t &_inverseproblem,
		const double _Delta,
		const unsigned int _maxiter,
		const unsigned int _maxinneriter,
		Database &_database,
		const enum DetermineStepWidthFactory::stepwidth_enumeration _stepwidth_type,
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
	C(0.9),
	stepwidth_type(_stepwidth_type)
{}

void LandweberMinimizer::setC(const double _C)
{
	if ((_C <= 0.)) // || ( _C > 1.))
		throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("C");

	const_cast<double&>(C) = _C;
}

GeneralMinimizer::ReturnValues
LandweberMinimizer::operator()(
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
//	const Norm & DualNormX = *DualSpaceX.getNorm();
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
	ReturnValues returnvalues;
	returnvalues.NumberOuterIterations = 0;
	// set iterate 'x' as start vector 'x0'
	returnvalues.m_solution = _problem->x->getSpace()->createElement();
	*returnvalues.m_solution = _startvalue;
	returnvalues.m_residual = _problem->y->getSpace()->createElement();
	// calculate starting residual and norm
	returnvalues.residuum = calculateResidual(
			_problem,
			returnvalues.m_residual);
	const double ynorm = NormY(y);

	// build data tuple for iteration information
	Table& per_iteration_table = database.addTable("per_iteration");
	Table::Tuple_t per_iteration_tuple = preparePerIterationTuple(
			NormX.getPvalue(), NormY.getPvalue(), 1, SpaceX.getDimension());

	// build data tuple for overall information
	Table& overall_table = database.addTable("overall");
	Table::Tuple_t overall_tuple = prepareOverallTuple(
			NormX.getPvalue(), NormY.getPvalue(), 1, SpaceX.getDimension());
	overall_tuple.insert( std::make_pair("runtime_matrix_vector_products", MatrixVectorProduct.getTiming()), Table::Data );
	overall_tuple.insert( std::make_pair("runtime_vector_vector_products", ScalarVectorProduct.getTiming()), Table::Data );

	/// -# check stopping criterion
	bool StopCriterion = false;
	StopCriterion = (fabs(returnvalues.residuum/ynorm) <= TolY);

	// calculate some values prior to loop
	SpaceElement_ptr_t dual_solution = DualSpaceX.createElement();
	*dual_solution = _dualstartvalue;
	BregmanDistance Delta_p(
			NormX,
			dynamic_cast<const PowerTypeDualityMapping &>(J_p),
			J_p.getPower(),
			ScalarVectorProduct);
	double old_distance = 0.;
	if (!_truesolution->isZero()) {
		old_distance = Delta_p(
				returnvalues.m_solution,
				_truesolution,
				dual_solution) + 1e4*BASSOTOLERANCE; // make sure its larger
		BOOST_LOG_TRIVIAL(debug)
				<< "Starting Bregman distance is " << old_distance;
	}

	// step width algorithm
	ResidualFunctional::calculateResidual_t residualizer =
			boost::bind(&LandweberMinimizer::calculateResidual,
					boost::cref(*this), _1, _2);
	DetermineStepWidth_ptr_t stepwidth =
			DetermineStepWidthFactory::createInstance(
					_problem,
					stepwidth_type,
					C,
					residualizer);

	/// -# loop over stopping criterion
	while (!StopCriterion) {
		per_iteration_tuple.replace( "iteration", (int)returnvalues.NumberOuterIterations);
		BOOST_LOG_TRIVIAL(debug)
				<< "#" << returnvalues.NumberOuterIterations
				<< " with residual of " << returnvalues.residuum;
		BOOST_LOG_TRIVIAL(debug)
			<< "#" << returnvalues.NumberOuterIterations << ": "
			<< "||Ax_n-y||/||y|| is " << returnvalues.residuum/ynorm;
		per_iteration_tuple.replace( "relative_residual", returnvalues.residuum/ynorm);
		// check that distance truely decreases
		if (!_truesolution->isZero()) {
			const double new_distance =
					Delta_p(
							returnvalues.m_solution,
							_truesolution,
							dual_solution);
			BOOST_LOG_TRIVIAL(debug)
				<< "#" << returnvalues.NumberOuterIterations << ": "
				<< "Delta_p^{x^*_n}(x_n,x) is "
				<< new_distance;
			per_iteration_tuple.replace( "bregman_distance", new_distance);
//			assert( old_distance > new_distance );
			old_distance = new_distance;
			const double new_error = NormX(returnvalues.m_solution-_truesolution);
			BOOST_LOG_TRIVIAL(debug)
				<< "#" << returnvalues.NumberOuterIterations << ": "
				<< "||x_n-x|| is " << new_error;
			per_iteration_tuple.replace( "error", new_error);
		}
		BOOST_LOG_TRIVIAL(trace)
				<< "x_n is " << returnvalues.m_solution;
		BOOST_LOG_TRIVIAL(trace)
				<< "R_n is " << returnvalues.m_residual;
		const SpaceElement_ptr_t Jw = j_r( returnvalues.m_residual );
		BOOST_LOG_TRIVIAL(trace)
				<< "j_r (residual) is " << Jw;
		const SpaceElement_ptr_t u = A_t * Jw;
		BOOST_LOG_TRIVIAL(trace)
				<< "u is " << u;

		// use step width used in theoretical proof
		// (F. Schöpfer, 11.4.2014) too conservative! Line search instead
		double alpha = (*stepwidth)(
				dual_solution,
				u,
				returnvalues.m_solution,
				returnvalues.m_residual,
				returnvalues.residuum,
				TolX,
				0.
				);
		per_iteration_tuple.replace( "stepwidth", alpha);
		BOOST_LOG_TRIVIAL(trace)
			<< "Step width is " << alpha;

		// iterate: J_p (x_{n+1})
		*dual_solution -= alpha * u;
		BOOST_LOG_TRIVIAL(trace)
				<< "x^*_n+1 is " << dual_solution;

		// finally map back from X^{\conj} to X: x_{n+1}
		*returnvalues.m_solution = J_q(dual_solution);
		BOOST_LOG_TRIVIAL(trace)
				<< "x_n+1 is " << returnvalues.m_solution;
		*_problem->x = returnvalues.m_solution;

		// update residual
		returnvalues.residuum = calculateResidual(
				_problem,
				returnvalues.m_residual);

		// check iterations count/wall time
		boost::chrono::high_resolution_clock::time_point timing_intermediate =
				boost::chrono::high_resolution_clock::now();
		++returnvalues.NumberOuterIterations;
		const double current_relative_residuum = fabs(returnvalues.residuum/ynorm);
		StopCriterion =
				CheckIterations(returnvalues.NumberOuterIterations)
				|| CheckResiduum(current_relative_residuum)
				|| CheckWalltime(boost::chrono::duration<double>(timing_intermediate - timing_start));

		// print intermediat solution
		printIntermediateSolution(
				 returnvalues.m_solution->getVectorRepresentation(),
				 A.getMatrixRepresentation(),
				returnvalues.NumberOuterIterations);

		// submit current tuple
		per_iteration_table.addTuple(per_iteration_tuple);
	}

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	std::cout << "The operation took " << boost::chrono::duration<double>(timing_end - timing_start) << "." << std::endl;

	// submit overall_tuple
	overall_tuple.replace( "iterations", returnvalues.NumberOuterIterations );
	overall_tuple.replace( "relative_residual", returnvalues.residuum );
	overall_tuple.replace( "runtime",
			boost::chrono::duration_cast<boost::chrono::duration<double> >(timing_end - timing_start).count() );
	overall_tuple.replace( "matrix_vector_products", (int)MatrixVectorProduct.getCount() );
	overall_tuple.replace( "vector_vector_products", (int)ScalarVectorProduct.getCount() );
	overall_tuple.replace( "runtime_matrix_vector_products", MatrixVectorProduct.getTiming() );
	overall_tuple.replace( "runtime_vector_vector_products", ScalarVectorProduct.getTiming() );
	overall_table.addTuple(overall_tuple);

	return returnvalues;
}
