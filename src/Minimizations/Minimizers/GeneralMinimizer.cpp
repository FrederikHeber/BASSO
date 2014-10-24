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

#include "Minimizations/Norms/NormFactory.hpp"

GeneralMinimizer::GeneralMinimizer(
		const DualityMappingsContainer &_container,
		const double _NormY,
		const double _PowerY,
		const double _Delta,
		const unsigned int _maxiter,
		Database &_database,
		const unsigned int _outputsteps
		) :
	DualityMappingsContainer(_container),
	val_NormY(_NormY),
	PowerY(_PowerY),
	Delta(_Delta),
	MaxOuterIterations(_maxiter),
	TolX(1e-6),
	TolY(Delta),
	TolFun(1e-12),
	outputsteps(_outputsteps),
	internal_NormX(NormFactory::createInstance(val_NormX)),
	internal_NormY(NormFactory::createInstance(val_NormY)),
	internal_DualNormX(NormFactory::createInstance(val_DualNormX)),
	NormX(*internal_NormX),
	NormY(*internal_NormY),
	DualNormX(*internal_DualNormX),
	j_r(val_NormY),
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
		<< "p is " << val_NormX
		<< ", q is " << val_DualNormX
		<< ", r is " << val_NormY
		<< ", power of J_p is " <<  PowerX
		<< ", power of J_q is " <<  DualPowerX
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
