/*
 * CommandLineOptions.hpp
 *
 *  Created on: Jan 27, 2015
 *      Author: heber
 */

#ifndef COMMANDLINEOPTIONS_HPP_
#define COMMANDLINEOPTIONS_HPP_

#include "BassoConfig.h"

#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
#include <string>
#include <vector>

#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "Minimizations/Minimizers/Searchspace/LastNSearchDirections.hpp"

/** This class contains the default command-line options that all of algorithms need.
 * It may be specialized to add specific commands for a specific executable.
 *
 * Under the hood we use boost::program_options to parse the options.
 */
class CommandLineOptions
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

	/** Shows the help conditions if desired.
	 *
	 * \param char array with program name, i.e. argv[0]
	 * \return true - help shown, so exit, false - else
	 */
	bool showHelpConditions(const char * const program_name) const;

	/** Sets the verbosity as specified by the command-line option.
	 *
	 */
	void setVerbosity() const;

	/** Checks whether the values make any sense.
	 *
	 * \return true - values sensible, false - not
	 */
	bool checkSensibility() const;

	/** Sets other values that require checkSensibility().
	 *
	 */
	void setSecondaryValues();

	//!> enumeration of all known dualities (Don't forget to add string literal to TypeNames)
	enum DualityContainerType {
		defaulttype=0,          //!< defaulttype
		regularizedl1norm=1,    //!< regularizedl1norm
		MAX_DualityContainerType//!< MAX_DualityContainerType
	};

protected:
	/** Override this function to add more options to \a desc.
	 *
	 */
	virtual void internal_init() {}

	/** Override this function to add parsing of your additional options.
	 *
	 */
	virtual void internal_parse() {}

	/** Override this function to add more conditions when help is shown.
	 *
	 */
	virtual bool internal_help_conditions() const { return false; }

	/** Override this function to add more sensibility checks.
	 *
	 */
	virtual bool internal_checkSensibility() const { return true; }

	/** Override this function to add setting secondary values after parsing.
	 *
	 */
	virtual void internal_setSecondaryValues() const {}

private:
	/** Internal function with own help condition checking.
	 *
	 * @return true - show help, false - else
	 */
	bool furtherHelpConditions() const;

protected:
	//!> container for all options
	boost::program_options::options_description desc;
	//!> key-value map of all parsed options
	boost::program_options::variables_map vm;

public:
	// primary options: set by command-line
	std::string algorithm_name;
	double C;
	bool calculateAngles;
	bool database_replace;
	double delta;
	bool enforceRandomMapping;
	bool inexactLinesearch;
	boost::filesystem::path iteration_file;
	unsigned int maxinneriter;
	std::string minlib;
	double normx;
	double normy;
	unsigned int N;
	unsigned int outputsteps;
	double powerx;
	double powery;
	double regularization_parameter;
	std::string searchspace_type;
	unsigned int stepwidth_type;
	double tau;
	std::vector<std::string> tuple_parameters;
	enum LastNSearchDirections::UpdateAlgorithmType updatetype;
	std::vector<double> wolfe_constants;

	// secondary options: set by other options
	DualityContainerType dualitytype;
	enum MinimizerFactory::InstanceType type;
};



#endif /* COMMANDLINEOPTIONS_HPP_ */
