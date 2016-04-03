/*
 * SolverFactory.hpp
 *
 *  Created on: Jan 28, 2015
 *      Author: heber
 */

#ifndef SOLVERFACTORY_SOLVERFACTORY_HPP_
#define SOLVERFACTORY_SOLVERFACTORY_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "Minimizations/types.hpp"

class CommandLineOptions;

/** This class contains static functions to ease in building executables
 * that require solving of inverse problems.
 *
 * The difficulty is that all classes in the folder Minimizations/
 * need many options but the Options/ folder should not be
 * known as the manner of parsing these options has nothing to do with
 * the inverse problem solving. Hence, we basically need an adaptor
 * between the CommandLineOptions and derived classes, wherein all the
 * required information for setting up and solving the the inverse
 * problem is contained, and the Minimization stuff, wherein all
 * functionality is contained for actually solving the problems.
 */
struct SolverFactory
{
	/** Creates the inverse problem given the information in \a _opts.
	 *
	 * \param _opts CommandLineOptions struct containing all information
	 * \param  _matrix matrix of the inverse problem
	 * \param  _rhs right-hand side of the problem
	 * \return constructed inverse problem instance
	 */
	static InverseProblem_ptr_t createInverseProblem(
			const CommandLineOptions &_opts,
			const Eigen::MatrixXd &_matrix,
			const Eigen::VectorXd &_rhs);

	/** Creates the inverse problem given the information in \a _opts.
	 *
	 * \param _opts CommandLineOptions struct containing all information
	 * \param _matrix_first_factor first factor of matrix of the inverse problem
	 * \param _matrix_second_factor second factor of matrix of the inverse problem
	 * \param _rhs right-hand side of the problem
	 * \return constructed inverse problem instance
	 */
	static InverseProblem_ptr_t createInverseProblemFromFactors(
			const CommandLineOptions &_opts,
			const Eigen::MatrixXd &_matrix_first_factor,
			const Eigen::MatrixXd &_matrix_second_factor,
			const Eigen::VectorXd &_rhs);

	/** Creates a database instance from the information contained in \a _opts.
	 *
	 * @param _opts CommandLineOptions struct containing all information
	 * @return constructed Database instance
	 */
	static Database_ptr_t createDatabase(
			const CommandLineOptions &_opts);

	/** Creates a minimizer instance from the information contained in \a _opts.
	 *
	 * @param _opts CommandLineOptions struct containing all information
	 * @param _inverseproblem inverse problem to be solved/minimized
	 * @param _database database to write iteration data to
	 * @return NULL - something went wrong,
	 * 		   else - constructed minimizer instance, ready to go
	 *
	 */
	static MinimizerFactory::instance_ptr_t createMinimizer(
		const CommandLineOptions &_opts,
		InverseProblem_ptr_t &_inverseproblem,
		Database_ptr_t &_database
		);
};



#endif /* SOLVERFACTORY_SOLVERFACTORY_HPP_ */
