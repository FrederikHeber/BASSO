
#include "BassoConfig.h"

#include <Eigen/Dense>
#include <fstream>
#include <string>

#include "CommandLineOptions/MatrixFactorizerOptions.hpp"
#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"


int main(int argc, char **argv)
{
	/// some required parameters
	MatrixFactorizerOptions opts;
	opts.init();

	// parse options
	opts.parse(argc, argv);

	if (opts.showHelpConditions(argv[0]))
		return 1;

	// set verbosity level
	opts.setVerbosity();

	if (!opts.checkSensibility())
		return 255;
	opts.setSecondaryValues();

	/// parse the matrices
	Eigen::MatrixXd data;
	Eigen::MatrixXd spectral_matrix;
	Eigen::MatrixXd picture_matrix;
	{
		using namespace MatrixIO;

		{
			std::ifstream ist(opts.data_file.string().c_str());
			if (ist.good())
				try {
					ist >> data;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully parse data matrix from " << opts.data_file.string() << std::endl;
					return 255;
				}
			else {
				std::cerr << "Failed to open " << opts.data_file.string() << std::endl;
				return 255;
			}

		}
	}
	// print parsed matrix and vector if small or high verbosity requested
	if ((data.innerSize() > 10) || (data.outerSize() > 10)) {
		BOOST_LOG_TRIVIAL(trace)
			<< "We solve for Y=K*X with Y = "
			<< data << "." << std::endl;
	} else {
		BOOST_LOG_TRIVIAL(info)
					<< "We solve for Y=K*X with Y = "
					<< data << "." << std::endl;
	}

	/// loop over pixel and channel dimensions

	// construct and solve (approximately) inverse problem

	/// output solution

	/// exit
	return 0;
}
