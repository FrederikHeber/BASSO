/*
 * AuxiliaryConstraintsFactory.cpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "AuxiliaryConstraintsFactory.hpp"

// boost::log uses spirit v3 and is incompatible with v2
#define BOOST_SPIRIT_USE_PHOENIX_V3 1
#include <boost/spirit/include/lex_lexertl.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>

#include <cassert>
#include <string>

#include "Log/Logging.hpp"

#include "Solvers/AuxiliaryConstraints/NonnegativeConstraint.hpp"
#include "Solvers/AuxiliaryConstraints/NonpositiveConstraint.hpp"
#include "Solvers/AuxiliaryConstraints/UnityConstraint.hpp"
#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraints_AND.hpp"

namespace lex = boost::spirit::lex;

enum token_ids
{
    ID_WORD = 1000,
    ID_EOL,
    ID_CHAR
};

template <typename Lexer>
struct boolean_combination_tokens : lex::lexer<Lexer>
{
	boolean_combination_tokens()
    {
        // define tokens (the regular expression to match and the corresponding
        // token id) and add them to the lexer
        this->self.add
            ("[^ \t\n]+", ID_WORD) // words (anything except ' ', '\t' or '\n')
            ("\n", ID_EOL)         // newline characters
            (".", ID_CHAR)         // anything else is a plain character
        ;
    }
};

struct parser
{
    typedef bool result_type;

    // the function operator gets called for each of the matched tokens
    // c, l, w are references to the counters used to keep track of the numbers
    template <typename Token>
    bool operator()(
    		Token const& t,
			AuxiliaryConstraintsFactory &f,
			AuxiliaryConstraints::ptr_t &criterion,
			AuxiliaryConstraintsFactory::Combination_types &combiner
			) const
    {
        switch (t.id()) {
        case ID_EOL:        // matched a newline character
        case ID_WORD:       // matched a word
        	// instaniate
        {
        	std::stringstream tokenword;
        	tokenword << t.value();
        	AuxiliaryConstraintsFactory::AuxiliaryConstraints_types stop_type =
        			f.getCriterionTypeForName(tokenword.str());
        	AuxiliaryConstraintsFactory::Combination_types combination_type =
        			f.getCombinationTypeForName(tokenword.str());
        	if (stop_type != AuxiliaryConstraintsFactory::MAX_AuxiliaryConstraints_types) {
            	BOOST_LOG_TRIVIAL(debug)
    				<< "Matched instance: " << tokenword.str();
        		AuxiliaryConstraints::ptr_t newcriterion =
        				f.createCriterion( stop_type );
        		if (combiner != AuxiliaryConstraintsFactory::NoCombination) {
        			criterion =
        					f.createCombination(combiner, criterion, newcriterion);
        			combiner = AuxiliaryConstraintsFactory::NoCombination;
        		} else {
        			criterion = newcriterion;
        		}
        	} else if (combination_type != AuxiliaryConstraintsFactory::MAX_Combination_types) {
            	BOOST_LOG_TRIVIAL(debug)
    				<< "Matched combination: " << tokenword.str();
        		// check whether it's an operator
        		assert( combiner == AuxiliaryConstraintsFactory::NoCombination);
        		combiner = combination_type;
        	} else {
        		BOOST_LOG_TRIVIAL(error)
        				<< "Cannot match instance: " << t.value();
        	}
            break;
        }
        case ID_CHAR:       // matched something else
//        	BOOST_LOG_TRIVIAL(info)
//				<< "Matched a char." << t.value();
        	// check whether it's a boolean operator
            break;
        }
        return true;        // always continue to tokenize
    }
};

AuxiliaryConstraintsFactory::AuxiliaryConstraintsFactory()
{
	// fill string to criterion type map
	StringToCriterionTypeMap["Nonnegative"] = Nonnegative;
	StringToCriterionTypeMap["Nonpositive"] = Nonpositive;
	StringToCriterionTypeMap["Unity"] = Unity;

	// fill string to combiner type map
	StringToCombinationTypeMap["&&"] = Combination_AND;
}

AuxiliaryConstraints::ptr_t
AuxiliaryConstraintsFactory::create(
		const std::string &_criteria_line)
{
	AuxiliaryConstraints::ptr_t criterion;
	Combination_types combiner = NoCombination;

    // these variables are used to count characters, words and lines
    // create the token definition instance needed to invoke the lexical analyzer
	boolean_combination_tokens<lex::lexertl::lexer<> > boolean_combination_functor;

    // tokenize the given string, the bound functor gets invoked for each of
    // the matched tokens
    char const* first = _criteria_line.c_str();
    char const* last = &first[_criteria_line.size()];
    bool r = lex::tokenize(first, last, boolean_combination_functor,
        boost::bind(parser(), _1,
        		boost::ref(*this),
				boost::ref(criterion),
				boost::ref(combiner)));

    if (!r) {
    	BOOST_LOG_TRIVIAL(error)
    			<< "Could not parse criteria line: " << _criteria_line;
    	criterion.reset();
    }

	return criterion;
}

AuxiliaryConstraints::ptr_t
AuxiliaryConstraintsFactory::createCriterion(
		const enum AuxiliaryConstraints_types _type)
{
	AuxiliaryConstraints::ptr_t criterion;

	switch (_type) {
	case Nonnegative:
		criterion.reset( new NonnegativeConstraint() );
		break;
	case Nonpositive:
		criterion.reset( new NonpositiveConstraint() );
		break;
	case Unity:
		criterion.reset( new UnityConstraint() );
		break;
	case MAX_AuxiliaryConstraints_types:
		criterion.reset();
		break;
	default:
		assert(0); // unknown criteria type
		break;
	}

	return criterion;
}

AuxiliaryConstraints::ptr_t
AuxiliaryConstraintsFactory::createCombination(
		const enum Combination_types _type,
		const AuxiliaryConstraints::ptr_t &_left,
		const AuxiliaryConstraints::ptr_t &_right)
{
	AuxiliaryConstraints::ptr_t criterion;

	switch (_type) {
	case Combination_AND:
		criterion.reset( new AuxiliaryConstraints_AND(_left, _right) );
		break;
	case NoCombination:
	case MAX_Combination_types:
		criterion.reset();
		break;
	default:
		assert(0); // unknown criteria type
		break;
	}

	return criterion;
}

enum AuxiliaryConstraintsFactory::AuxiliaryConstraints_types
AuxiliaryConstraintsFactory::getCriterionTypeForName(const std::string &_name)
{
	StringToCriterionTypeMap_t::const_iterator iter = StringToCriterionTypeMap.find(_name);
	if (iter == StringToCriterionTypeMap.end())
		return MAX_AuxiliaryConstraints_types;
	else
		return iter->second;
}

enum AuxiliaryConstraintsFactory::Combination_types
AuxiliaryConstraintsFactory::getCombinationTypeForName(const std::string &_name)
{
	StringToCombinationTypeMap_t::const_iterator iter = StringToCombinationTypeMap.find(_name);
	if (iter == StringToCombinationTypeMap.end())
		return MAX_Combination_types;
	else
		return iter->second;
}
