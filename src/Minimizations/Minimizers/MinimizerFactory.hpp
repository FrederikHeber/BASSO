/*
 * MinimizerFactory.hpp
 *
 *  Created on: Apr 28, 2014
 *      Author: heber
 */

#ifndef MINIMIZERFACTORY_HPP_
#define MINIMIZERFACTORY_HPP_

#include "BassoConfig.h"

#include <boost/shared_ptr.hpp>
#include <string>

#include "Minimizations/types.hpp"
#include "Minimizations/Minimizers/Searchspace/LastNSearchDirections.hpp"
#include "Minimizations/Minimizers/StepWidths/DetermineStepWidthFactory.hpp"

class CommandLineOptions;
class Database;
class DualityMappingsContainer;
class GeneralMinimizer;

/** This class is a factory for GeneralMinimizers, producing
 * specific instances but handing back only the general interface
 * such that they can be used in the same way.
 *
 */
class MinimizerFactory
{
public:
	typedef boost::shared_ptr<GeneralMinimizer> instance_ptr_t;

	MinimizerFactory() {}
	~MinimizerFactory() {}

	//!> enumeration of all known instances (Don't forget to add string literal to TypeNames)
	enum InstanceType {
		landweber=0,
		sequentialsubspace=1,
		sequentialsubspace_noise=2,
		MAX_InstanceType
	};

	/** Produces the desired instance.
	 *
	 * @param _opts CommandLineOptions struct containing all information
	 * @param _inverseproblem inverse problem to solve
	 * @param _database database to store iteration information to
	 * @return wrapped instance of desired \a _type
	 */
	instance_ptr_t createInstance(
			const CommandLineOptions &_opts,
			const InverseProblem_ptr_t &_inverseproblem,
			Database &_database
			);

	/** Helper function to check whether the given \a _name designates a
	 * known type.
	 *
	 * @param _name name of InstanceType to check
	 * \return true - name ok, false - name unknown
	 */
	static bool isValidTypeName(
			const std::string &_name)
	{
		return getTypeNamesIndex(_name) != MAX_InstanceType;
	}


	/** Helper function to resolve internal type id for a given type name.
	 *
	 * @param _name name of instance type
	 * @return id of this instance type
	 */
	static enum InstanceType getTypeForName(
			const std::string &_name);

	/** Helper function to get name to an internal instance type.
	 *
	 * In contrast to direct access to TypeNames we also check that
	 * _type is valid.
	 *
	 * @param _type id of instance type
	 * @return string containing type's name
	 */
	static const std::string& getNameForType(
			const enum InstanceType &_type);

protected:
	/** Internal function to look up the index for the given name.
	 *
	 * @param _name name of InstanceType to look up in TypeNames
	 * @return index in array
	 */
	static unsigned int getTypeNamesIndex(
			const std::string &_name);

public:
	//!> name of each known instance type (Don't forget to add enum to InstanceType)
	static const std::string TypeNames[];
};


#endif /* MINIMIZERFACTORY_HPP_ */
