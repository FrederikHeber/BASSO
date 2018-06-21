#include <Eigen/Dense>

static Eigen::VectorXd getRandomVector(
                const int _N)
{
        Eigen::VectorXd random_vector(_N);
        random_vector.setRandom();
        return random_vector;
}

Eigen::MatrixXd getRandomMatrix(
                const int _N, /* target dim */
                const int _M /* source dim */
                )
{
        Eigen::MatrixXd matrix(_N,_M);
        matrix.setRandom();
        return matrix;
}

