/*
 * CommandLineOptions.hpp
 *
 *  Created on: Jan 27, 2015
 *      Author: heber
 */

#ifndef COMMANDLINEOPTIONS_HPP_
#define COMMANDLINEOPTIONS_HPP_

#include "BassoConfig.h"

#include <iosfwd>
#include <string>
#include <vector>

#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "Minimizations/Minimizers/Searchspace/LastNSearchDirections.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingArguments.hpp"
#include "Options/Options.hpp"
#include "Options/filesystem_path_serialization.hpp"

/** This class contains the default command-line options that all of algorithms need.
 * It may be specialized to add specific commands for a specific executable.
 *
 * Under the hood we use boost::program_options to parse the options.
 */
class CommandLineOptions : public Options
{
public:
	CommandLineOptions();
	virtual ~CommandLineOptions();

	/** Initializes the options.
	 *
	 */
	void init();

	/** Performs the actual parsing.
	 *
	 * \param argc argument count
	 * \param **argv array of arguments
	 */
	void parse(int argc, char **argv);

	/** Checks whether the values make any sense.
	 *
	 * \return true - values sensible, false - not
	 */
	bool checkSensibility() const;

	/** Sets other values that require checkSensibility().
	 *
	 */
	void setSecondaryValues();

	/** Stores all variables as "key = value" pairs per line under the given
	 * stream \a _output.
	 *
	 * @param _output output stream to write to
	 */
	void store(std::ostream &_output) const;

protected:
	/** Override this function to add more options to \a desc.
	 *
	 */
	virtual void internal_init() {}

	/** Override this function to add parsing of your additional options.
	 *
	 */
	virtual void internal_parse() {}

	/** Override this function to add more sensibility checks.
	 *
	 */
	virtual bool internal_checkSensibility() const { return true; }

	/** Override this function to add setting secondary values after parsing.
	 *
	 */
	virtual void internal_setSecondaryValues() const {}

	/** Override this function to store values from deriving classes.
	 *
	 */
	virtual void internal_store(std::ostream &_output) const {}

private:
	bool checkSensibility_config() const;
	bool checkSensibility_delta() const;
	bool checkSensibility_everynthtuple() const;
	bool checkSensibility_OrthogonalDirections() const;
	bool checkSensibility_regularizationparameter() const;
	bool checkSensibility_tau() const;
	bool checkSensibility_tuple_parameters() const;
	bool checkSensibility_algorithm() const;
	bool checkSensibility_minlib() const;
	bool checkSensibility_norms() const;
	bool checkSensibility_pvalues() const;
	bool checkSensibility_searchspace() const;
	bool checkSensibility_stepwidth_algorithm() const;
	bool checkSensibility_stopping_criteria() const;
	bool checkSensibility_updatealgorithm() const;
	bool checkSensibility_wolfeconstants() const;

public:
	// primary options: set by command-line
	std::string algorithm_name;
	std::string auxiliary_constraints;
	double C;
	bool calculateAngles;
	bool database_replace;
	double delta;
	bool enforceRandomMapping;
	unsigned int everynthtuple;
	bool inexactLinesearch;
	boost::filesystem::path iteration_file;
	unsigned int maxinneriter;
	unsigned int maxiter;
	unsigned int max_sfp_loops;
	double maxwalltime;
	std::string minlib;
	std::string type_spacex;
	std::string type_spacey;
	double px;
	double py;
	unsigned int N;
	enum LastNSearchDirections::OrthogonalizationType orthogonalization_type;
	unsigned int outputsteps;
	double powerx;
	double powery;
	double regularization_parameter;
	std::string searchspace_type;
	unsigned int stepwidth_type;
	std::string stopping_criteria;
	double tau;
	double tolerance_linesearch;
	double tolerance_spacex;
	std::vector<std::string> tuple_parameters;
	enum LastNSearchDirections::UpdateAlgorithmType updatetype;
	unsigned int verbose;
	std::vector<double> wolfe_constants;

	// secondary options: set by other options
	enum MinimizerFactory::InstanceType type;

private:
	friend class boost::serialization::access;
	template<class Archive>
	void save(Archive & ar, const unsigned int version) const
	{
		// CommandLineOptions
		ar & boost::serialization::base_object<const Options>(*this);

		ar & algorithm_name;
		ar & auxiliary_constraints;
		ar & C;
		ar & calculateAngles;
		ar & database_replace;
		ar & delta;
		ar & enforceRandomMapping;
		ar & everynthtuple;
		ar & inexactLinesearch;
		ar & iteration_file;
		ar & maxinneriter;
		ar & maxiter;
		ar & max_sfp_loops;
		ar & maxwalltime;
		ar & minlib;
		ar & type_spacex;
		ar & type_spacey;
		ar & px;
		ar & py;
		ar & N;
		ar & orthogonalization_type;
		ar & outputsteps;
		ar & powerx;
		ar & powery;
		ar & regularization_parameter;
		ar & searchspace_type;
		ar & stepwidth_type;
		ar & stopping_criteria;
		ar & tau;
		ar & tolerance_linesearch;
		ar & tolerance_spacex;
		ar & tuple_parameters;
		ar & updatetype;
		ar & wolfe_constants;
	}

	template<class Archive>
	void load(Archive & ar, const unsigned int version)
	{
		// CommandLineOptions
		ar & boost::serialization::base_object<Options>(*this);

		ar & algorithm_name;
		ar & auxiliary_constraints;
		ar & C;
		ar & calculateAngles;
		ar & database_replace;
		ar & delta;
		ar & enforceRandomMapping;
		ar & everynthtuple;
		ar & inexactLinesearch;
		ar & iteration_file;
		ar & maxinneriter;
		ar & maxiter;
		ar & max_sfp_loops;
		ar & maxwalltime;
		ar & minlib;
		ar & type_spacex;
		ar & type_spacey;
		ar & px;
		ar & py;
		ar & N;
		ar & orthogonalization_type;
		ar & outputsteps;
		ar & powerx;
		ar & powery;
		ar & regularization_parameter;
		ar & searchspace_type;
		ar & stepwidth_type;
		ar & stopping_criteria;
		ar & tau;
		ar & tolerance_linesearch;
		ar & tolerance_spacex;
		ar & tuple_parameters;
		ar & updatetype;
		ar & wolfe_constants;

		// set secondary values
		setSecondaryValues();
	}

	BOOST_SERIALIZATION_SPLIT_MEMBER()
};



#endif /* COMMANDLINEOPTIONS_HPP_ */
