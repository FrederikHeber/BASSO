#include "random_matrix_inverseproblem_1.hpp"

int main(int argc, char **argv) {
	const int innerSize = 500;
	const int outerSize = 100;
	const Eigen::MatrixXd matrix = getRandomMatrix(outerSize,innerSize);
	const Eigen::VectorXd rhs = getRandomVector(outerSize);

	// place value of p and value of power type into boost any list
	InverseProblemFactory::args_t args_SpaceX;
	InverseProblemFactory::args_t args_SpaceY; 
	args_SpaceX += boost::any(1.1), boost::any(2.);
	args_SpaceY += boost::any(2.), boost::any(2.);

	// first, create the spaces X and Y
	NormedSpace_ptr_t X = NormedSpaceFactory::create(
		        innerSize, "lp", args_SpaceX);
	NormedSpace_ptr_t Y = NormedSpaceFactory::create(
		        outerSize, "lp", args_SpaceY);

	// then create the right-hand side in the space Y
	SpaceElement_ptr_t y = ElementCreator::create(Y, rhs);

	// and the LinearMapping A: X -> Y
	Mapping_ptr_t A = LinearMappingFactory::createInstance(X,Y, matrix, false);

	// finally the inverse problem Ax=y
	const InverseProblem_ptr_t inverseproblem(new InverseProblem(A,X,Y,y) );

	// create options for the solver
	CommandLineOptions opts;
	opts.delta = 1e-4; /* default value */
	opts.maxiter = 50; /* default value */
	opts.orthogonalization_type = LastNSearchDirections::MetricOrthogonalization;
	// stop when max iterations are exceeded or relative residum change
	// is less than threshold
	opts.stopping_criteria = "MaxIterationCount || RelativeResiduum";
	opts.setVerbosity();
	opts.setSecondaryValues();

	// create a dummy database (which takes all iteration-related information
	// such as timings, number of various function calls, ...) but does not store
	// them unless we use setDatabaseFile()
	Database_ptr_t database(new SQLDatabase());

	// create minimizer/solver SESOP with one direction
	SequentialSubspaceMinimizer *solver = new SequentialSubspaceMinimizer(opts,
		inverseproblem,
		*database);
	solver->setN(1);

	// create starting point as zero in X
	SpaceElement_ptr_t solution_start = inverseproblem->x->getSpace()->createElement();
	SpaceElement_ptr_t zero = inverseproblem->x->getSpace()->createElement();

	// create dual starting point as zero in X^\ast
	SpaceElement_ptr_t dual_start = inverseproblem->x->getSpace()->getDualSpace()->createElement();

	// solve Ax=y using SESOP
	GeneralMinimizer::ReturnValues result = (*solver)( 
		inverseproblem,
		solution_start,
		dual_start,
		zero );

	// check whether solving was successful
	if (result.status == GeneralMinimizer::ReturnValues::error) {
		LOG(error, "Something went wrong during minimization.");
	}

	// print iterations and remaining residual
	std::cout << "Solution after " << result.NumberOuterIterations
		<< " with relative residual of " << result.residuum/inverseproblem->y->Norm() 
		<< std::endl;
};