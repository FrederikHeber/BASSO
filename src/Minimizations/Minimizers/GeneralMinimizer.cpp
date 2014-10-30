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

#include <boost/bind.hpp>
#include <boost/log/trivial.hpp>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Norms/NormFactory.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMappingFactory.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

GeneralMinimizer::GeneralMinimizer(
		const InverseProblem_ptr_t &_inverseproblem,
		const double _Delta,
		const unsigned int _maxiter,
		Database &_database,
		const unsigned int _outputsteps
		) :
	val_NormY(_inverseproblem->y->getSpace()->getNorm()->getPvalue()),
	PowerY(_inverseproblem->y->getSpace()->getDualityMapping()->getPower()),
	Delta(_Delta),
	MaxOuterIterations(_maxiter),
	TolX(1e-6),
	TolY(Delta),
	TolFun(1e-12),
	outputsteps(_outputsteps),
	internal_NormX(_inverseproblem->x->getSpace()->getNorm()),
	internal_NormY(_inverseproblem->y->getSpace()->getNorm()),
	internal_DualNormX(_inverseproblem->x->getSpace()->getDualSpace()->getNorm()),
	internal_J_p(_inverseproblem->x->getSpace()->getDualityMapping()),
	internal_J_q(
			_inverseproblem->x->getSpace()->getDualSpace()->getDualityMapping()
	),
	internal_j_r(_inverseproblem->y->getSpace()->getDualityMapping()),
	NormX(*internal_NormX),
	NormY(*internal_NormY),
	DualNormX(*internal_DualNormX),
	J_p(*internal_J_p),
	J_q(*internal_J_q),
	j_r(*internal_j_r),
	database(_database),
	matrix_vector_fctor(
			boost::bind(
					static_cast<const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type
						(Eigen::MatrixBase<Eigen::MatrixXd>::*)(const Eigen::MatrixBase<Eigen::VectorXd>&) const>(
								&Eigen::MatrixBase<Eigen::MatrixXd>::operator*),
								_1, _2
			)
	),
	MatrixVectorProduct(matrix_vector_fctor),
	scalar_vector_fctor(
			boost::bind(
					static_cast<Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType
						(Eigen::MatrixBase<Eigen::VectorXd>::*)(const Eigen::MatrixBase<Eigen::VectorXd>&) const>(
								&Eigen::MatrixBase<Eigen::VectorXd>::dot),
								_1, _2
			)
	),
	ScalarVectorProduct(scalar_vector_fctor)
{
	BOOST_LOG_TRIVIAL(debug)
		<< "p is " << _inverseproblem->x->getSpace()->getNorm()->getPvalue()
		<< ", q is " << _inverseproblem->x->getSpace()->getDualSpace()->getNorm()->getPvalue()
		<< ", r is " << val_NormY
		<< ", power of J_p is " <<  _inverseproblem->x->getSpace()->getDualityMapping()->getPower()
		<< ", power of J_q is " <<  _inverseproblem->x->getSpace()->getDualSpace()->getDualityMapping()->getPower()
		<< ", power of J_r is " <<  PowerY;

	// set tolerances values
	J_p.setTolerance(TolX);
	J_q.setTolerance(TolX);
	j_r.setTolerance(TolY);
}

double GeneralMinimizer::calculateResidual(
		const Eigen::VectorXd &_x0,
		const Eigen::MatrixXd &_A,
		const Eigen::VectorXd &_y,
		Eigen::VectorXd &_residual
		) const
{
	_residual = MatrixVectorProduct(_A,_x0) - _y;
	return NormY(_residual);
}

double GeneralMinimizer::calculateResidual(
		const InverseProblem_ptr_t &_problem,
		SpaceElement_ptr_t &_residual
		) const
{
	*_residual =  MatrixVectorProduct(
			dynamic_cast<const LinearMapping &>(*_problem->A).getMatrixRepresentation(),
			_problem->x->getVectorRepresentation());
	*_residual -= _problem->y;
	return NormY(_residual);
}

void GeneralMinimizer::printIntermediateSolution(
		const Eigen::VectorXd &_solution,
		const Eigen::MatrixXd &_A,
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
