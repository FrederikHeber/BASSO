/*
 * StoppingCriteriaFactory.cpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "StoppingCriteriaFactory.hpp"

// boost::log uses spirit v3 and is incompatible with v2
#define BOOST_SPIRIT_USE_PHOENIX_V3 1
#include <boost/spirit/include/lex_lexertl.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>

#include <cassert>
#include <string>

#include "Log/Logging.hpp"

#include "Minimizations/Minimizers/StoppingCriteria/CheckIterationCount.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/CheckRelativeChangeResiduum.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/CheckRelativeResiduum.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/CheckResiduum.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/CheckWalltime.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingArguments.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_AND.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_NOT.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_OR.hpp"

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
			StoppingCriteriaFactory &f,
			StoppingCriterion::ptr_t &criterion,
			StoppingCriteriaFactory::Combination_types &combiner,
			const StoppingArguments &args
			) const
    {
        switch (t.id()) {
        case ID_EOL:        // matched a newline character
        case ID_WORD:       // matched a word
        	// instaniate
        {
        	std::stringstream tokenword;
        	tokenword << t.value();
        	StoppingCriteriaFactory::StoppingCriteria_types stop_type =
        			f.getCriterionTypeForName(tokenword.str());
        	StoppingCriteriaFactory::Combination_types combination_type =
        			f.getCombinationTypeForName(tokenword.str());
        	if (stop_type != StoppingCriteriaFactory::MAX_StoppingCriteria_types) {
            	BOOST_LOG_TRIVIAL(debug)
    				<< "Matched instance: " << tokenword.str();
        		StoppingCriterion::ptr_t newcriterion =
        				f.createCriterion( stop_type, args );
        		if (combiner != StoppingCriteriaFactory::NoCombination) {
        			criterion =
        					f.createCombination(combiner, criterion, newcriterion);
        			combiner = StoppingCriteriaFactory::NoCombination;
        		} else {
        			criterion = newcriterion;
        		}
        	} else if (combination_type != StoppingCriteriaFactory::MAX_Combination_types) {
            	BOOST_LOG_TRIVIAL(debug)
    				<< "Matched combination: " << tokenword.str();
        		// check whether it's an operator
        		assert( combiner == StoppingCriteriaFactory::NoCombination);
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

StoppingCriteriaFactory::StoppingCriteriaFactory()
{
	// fill string to criterion type map
	StringToCriterionTypeMap["MaxIterationCount"] = IterationCount;
	StringToCriterionTypeMap["RelativeChangeResiduum"] = RelativeChangeResiduum;
	StringToCriterionTypeMap["RelativeResiduum"] = RelativeResiduum;
	StringToCriterionTypeMap["Residuum"] = Residuum;
	StringToCriterionTypeMap["MaxWalltime"] = Walltime;

	// fill string to combiner type map
	StringToCombinationTypeMap["&&"] = Combination_AND;
	StringToCombinationTypeMap["||"] = Combination_OR;
	StringToCombinationTypeMap["!"] = Combination_NOT;
}

StoppingCriterion::ptr_t
StoppingCriteriaFactory::create(
		const std::string &_criteria_line,
		const StoppingArguments &_args)
{
	StoppingCriterion::ptr_t criterion;
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
				boost::ref(combiner),
				boost::cref(_args)));

    if (!r) {
    	BOOST_LOG_TRIVIAL(error)
    			<< "Could not parse criteria line: " << _criteria_line;
    	criterion.reset();
    }

	return criterion;
}

StoppingCriterion::ptr_t
StoppingCriteriaFactory::createCriterion(
		const enum StoppingCriteria_types _type,
		const StoppingArguments &_args)
{
	StoppingCriterion::ptr_t criterion;

	switch (_type) {
	case IterationCount:
		criterion.reset( new CheckIterationCount(_args) );
		break;
	case RelativeChangeResiduum:
		criterion.reset( new CheckRelativeChangeResiduum(_args) );
		break;
	case RelativeResiduum:
		criterion.reset( new CheckRelativeResiduum(_args) );
		break;
	case Residuum:
		criterion.reset( new CheckResiduum(_args) );
		break;
	case Walltime:
		criterion.reset( new CheckWalltime(_args) );
		break;
	case MAX_StoppingCriteria_types:
		criterion.reset();
		break;
	default:
		assert(0); // unknown criteria type
		break;
	}

	return criterion;
}

StoppingCriterion::ptr_t
StoppingCriteriaFactory::createCombination(
		const enum Combination_types _type,
		const StoppingCriterion::ptr_t &_left,
		const StoppingCriterion::ptr_t &_right)
{
	StoppingCriterion::ptr_t criterion;

	switch (_type) {
	case Combination_AND:
		criterion.reset( new StoppingCriterion_AND(_left, _right) );
		break;
	case Combination_NOT:
		// NOT operator comes before criterion, hence _right contains something
		criterion.reset( new StoppingCriterion_NOT(_right) );
		break;
	case Combination_OR:
		criterion.reset( new StoppingCriterion_OR(_left, _right) );
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

enum StoppingCriteriaFactory::StoppingCriteria_types
StoppingCriteriaFactory::getCriterionTypeForName(const std::string &_name)
{
	StringToCriterionTypeMap_t::const_iterator iter = StringToCriterionTypeMap.find(_name);
	if (iter == StringToCriterionTypeMap.end())
		return MAX_StoppingCriteria_types;
	else
		return iter->second;
}

enum StoppingCriteriaFactory::Combination_types
StoppingCriteriaFactory::getCombinationTypeForName(const std::string &_name)
{
	StringToCombinationTypeMap_t::const_iterator iter = StringToCombinationTypeMap.find(_name);
	if (iter == StringToCombinationTypeMap.end())
		return MAX_Combination_types;
	else
		return iter->second;
}
