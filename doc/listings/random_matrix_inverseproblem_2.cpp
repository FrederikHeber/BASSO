#include "random_matrix_inverseproblem_2.hpp"

int main(int argc, char **argv) {
	const int innerSize = 500;
	const int outerSize = 100;
	const Eigen::MatrixXd matrix =
		getRandomMatrix(outerSize,innerSize);
	const Eigen::VectorXd rhs = getRandomVector(outerSize);

	// create options for the solver
	CommandLineOptions opts;
	opts.algorithm_name = "SESOP";
	opts.delta = 1e-4; /* default value */
	opts.maxiter = 50; /* default value */
	opts.orthogonalization_type =
		LastNSearchDirections::MetricOrthogonalization;
	opts.px = 1.1;
	opts.py = 2.; /* already the default value */
	opts.powerx = 2.; /* already the default value */
	opts.powery = 2.; /* already the default value */
	// stop when max iterations are exceeded or relative residum
	// change is less than threshold
	opts.stopping_criteria = "MaxIterationCount || RelativeResiduum";
	opts.setVerbosity();
	opts.setSecondaryValues();

	// create a dummy database (which takes all iteration-related
	// information such as timings, number of various function
	// calls, ...) but does not store them unless we use
	// its setDatabaseFile()
	Database_ptr_t database =
		SolverFactory::createDatabase(opts);

	// prepare inverse problem
	InverseProblem_ptr_t inverseproblem =
		SolverFactory::createInverseProblem(
			opts, matrix, rhs);

	// solving
	InverseProblemSolver solver(
		inverseproblem,
		database,
		opts,
		false /* true solution calculation */);

	// create starting point as zero in X
	SpaceElement_ptr_t solution_start =
		inverseproblem->x->getSpace()->createElement();

	// solve Ax=y using SESOP
	GeneralMinimizer::ReturnValues result = solver(
		solution_start );

	// check whether solving was successful
	if (result.status == GeneralMinimizer::ReturnValues::error) {
		LOG(error, "Something went wrong during minimization.");
	}

	// print iterations and remaining residual
	std::cout << "Solution after " << result.NumberOuterIterations
		<< " with relative residual of "
		<< result.residuum/inverseproblem->y->Norm()
		<< std::endl;
};
