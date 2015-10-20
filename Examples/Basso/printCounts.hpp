/*
 * printCounts.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef BASSOBASE_PRINTCOUNTS_HPP_
#define BASSOBASE_PRINTCOUNTS_HPP_

#include "BassoConfig.h"

#include <vector>

#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/types.hpp"

struct printCounts
{
	printCounts(const std::vector<NormedSpace_ptr_t> &_spaces);

	template <class _type>
	void operator()(_type &) {
		int sum = 0;
		for (std::vector<NormedSpace_ptr_t>::const_iterator iter = spaces.begin();
				iter != spaces.end(); ++iter)
			sum += VectorSpaceOperations::getCountTiming<_type>((*iter)->getOpCounts().instance).first;
		std::cout << "\t" << VectorSpaceOperations::TypeToString<_type>()() << ": " << sum << std::endl;
	}
private:
	const std::vector<NormedSpace_ptr_t> &spaces;
};





#endif /* BASSOBASE_PRINTCOUNTS_HPP_ */
