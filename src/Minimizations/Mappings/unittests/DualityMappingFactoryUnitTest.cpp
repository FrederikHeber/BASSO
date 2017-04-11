/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*
 * DualityMappingFactoryUnitTest.cpp
 *
 *  Created on: Jul 01, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include "DualityMappingFactoryUnitTest.hpp"

#include <boost/bind.hpp>

#include "Minimizations/Mappings/DualityMappingFactory.hpp"
#include "Minimizations/Norms/NormFactory.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DualityMappingFactoryUnitTest );


void DualityMappingFactoryUnitTest::setUp()
{
}


void DualityMappingFactoryUnitTest::tearDown()
{
}

void DualityMappingFactoryUnitTest::equivalentTokensTest()
{
	// check whether the tokens in NormFactory and DualityMappingFactory
	// coincide
	typedef std::vector<std::string> tokens_t;

	const NormFactory::TokenCreatorMap_t NormFactory_TokenCreatorMap =
			NormFactory::getMap();
	tokens_t NormFactory_tokens(NormFactory_TokenCreatorMap.size());
	std::transform(
			NormFactory_TokenCreatorMap.begin(), NormFactory_TokenCreatorMap.end(),
			NormFactory_tokens.begin(),
			bind( &NormFactory::TokenCreatorMap_t::value_type::first, _1) );

	const DualityMappingFactory::TokenCreatorMap_t DualityMappingFactory_TokenCreatorMap =
			DualityMappingFactory::getMap();
	tokens_t DualityMappingFactory_tokens(DualityMappingFactory_TokenCreatorMap.size());;
	std::transform(
			DualityMappingFactory_TokenCreatorMap.begin(), DualityMappingFactory_TokenCreatorMap.end(),
			DualityMappingFactory_tokens.begin(),
			bind( &DualityMappingFactory::TokenCreatorMap_t::value_type::first, _1) );

	std::sort(NormFactory_tokens.begin(), NormFactory_tokens.end());
	std::sort(DualityMappingFactory_tokens.begin(), DualityMappingFactory_tokens.end());

	CPPUNIT_ASSERT( NormFactory_tokens == DualityMappingFactory_tokens );
}
