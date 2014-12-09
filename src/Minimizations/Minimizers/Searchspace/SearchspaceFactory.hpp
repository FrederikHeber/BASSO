/*
 * SearchspaceFactory.hpp
 *
 *  Created on: Nov 18, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_SEARCHSPACEFACTORY_HPP_
#define MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_SEARCHSPACEFACTORY_HPP_

#include "BassoConfig.h"

#include <boost/bimap.hpp>

#include "Minimizations/Minimizers/Searchspace/Searchspace.hpp"

/** This class produces the desired search space type.
 *
 */
class SearchspaceFactory
{
public:
	//!> enumeration of all known types
	enum SearchspaceType {
		InvalidType=0,
		LastNDirections,
		Nemirovsky,
		MAX_SearchspaceType
	};

	/** Creates a Searchspace object of the current \a InstanceType.
	 *
	 * @param _SearchDirectionSpace_ptr search direction space (for checks)
	 * @param _N number of search directions
	 * @param _MatrixVectorProduct_subspace counts for matrix-vector
	 * @param _ScalarVectorProduct_subspace counts for vector-vector
	 * @return shared_ptr containing created instance
	 */
	static Searchspace::ptr_t create(
			const NormedSpace_ptr_t &_SearchDirectionSpace_ptr,
			const unsigned int _N
			);

	/** Setter for the desired type to produce.
	 *
	 * @param _type type to produce
	 */
	static void setCurrentType(const enum SearchspaceType _type);

	//!> typedef for the map containing name and type associations
	typedef boost::bimap<enum SearchspaceType, std::string> NameTypeMap_t;

	static const std::string getName(const enum SearchspaceType _type);

	static bool isValidName(const std::string &_name);

	static const enum SearchspaceType getType(const std::string &_name);

	static bool isValidType(const enum SearchspaceType _type);

private:
	static void fillNameTypeMap();

private:
	//!> map containing name and type associations
	static NameTypeMap_t NameTypeMap;

	//!> current type to create
	static enum SearchspaceType InstanceType;
};



#endif /* MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_SEARCHSPACEFACTORY_HPP_ */
