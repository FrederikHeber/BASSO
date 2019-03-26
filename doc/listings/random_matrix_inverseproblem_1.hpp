#include <iostream>

#include <boost/assign.hpp>

#include <Eigen/Dense>

#include <Database/SQLDatabase.hpp>
#include <Log/Logging.hpp>
#include <Minimizations/Elements/ElementCreator.hpp>
#include <Minimizations/Elements/SpaceElement.hpp>
#include <Minimizations/InverseProblems/InverseProblem.hpp>
#include <Minimizations/InverseProblems/InverseProblemFactory.hpp>
#include <Minimizations/Mappings/LinearMapping.hpp>
#include <Minimizations/Mappings/MappingFactory.hpp>
#include <Minimizations/Minimizers/Searchspace/LastNSearchDirections.hpp>
#include <Minimizations/Minimizers/SequentialSubspaceMinimizer.hpp>
#include <Minimizations/Spaces/NormedSpaceFactory.hpp>
#include <Options/CommandLineOptions.hpp>
#include <Solvers/SolverFactory/SolverFactory.hpp>

#include "random_matrix.hpp"

using namespace boost::assign;

