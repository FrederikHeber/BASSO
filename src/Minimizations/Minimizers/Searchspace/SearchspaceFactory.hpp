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
};



#endif /* MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_SEARCHSPACEFACTORY_HPP_ */
