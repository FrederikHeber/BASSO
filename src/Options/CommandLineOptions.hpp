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
#include <iosfwd>
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

	/** Shows a help message in case an error parsing the option occurs.
	 *
	 */
	void showHelpinErrorCase() const;

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

	/** Stores all variables as "key = value" pairs per line under the given
	 * stream \a _output.
	 *
	 * @param _output output stream to write to
	 */
	void store(std::ostream &_output);

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

	/** Override this function to store values from deriving classes.
	 *
	 */
	virtual void internal_store(std::ostream &_output) const {}

private:
	bool checkSensibility_delta() const;
	bool checkSensibility_OrthogonalDirections() const;
	bool checkSensibility_regularizationparameter() const;
	bool checkSensibility_tau() const;
	bool checkSensibility_tuple_parameters() const;
	bool checkSensibility_algorithm() const;
	bool checkSensibility_minlib() const;
	bool checkSensibility_searchspace() const;
	bool checkSensibility_updatealgorithm() const;
	bool checkSensibility_wolfeconstants() const;


protected:
	/** Helper template to write "key = value" pair to a stream.
	 *
	 * Note that the template class \a T designates the type of the
	 * value.
	 *
	 * \note This function resides here to allow use by derived classes
	 * while not requiring to have a new namespace (i.e. we just use
	 * CommandLineOptions as a namespace).
	 *
	 * @param _output output stream
	 * @param _vm map containing all option keys and variables
	 * @param _token key as string
	 */
	template <class T>
	static
	void writeValue(
			std::ostream &_output,
			const boost::program_options::variables_map &_vm,
			const std::string &_token)
	{
		if (_vm.count(_token))
			_output << "\t" <<
				_token << " = " << _vm[_token].as<T>() << std::endl;
	}

protected:
	//!> container for all options combined
	boost::program_options::options_description desc_all;

	//!> container for options specifying the configuration file
	boost::program_options::options_description desc_config;

	//!> key-value map of all parsed options
	boost::program_options::variables_map vm;

public:
	// primary options: set by command-line
	std::string algorithm_name;
	double C;
	bool calculateAngles;
	boost::filesystem::path config_filename;
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
	enum LastNSearchDirections::OrthogonalizationType orthogonalization_type;
	unsigned int outputsteps;
	double powerx;
	double powery;
	double regularization_parameter;
	std::string searchspace_type;
	unsigned int stepwidth_type;
	double tau;
	std::vector<std::string> tuple_parameters;
	enum LastNSearchDirections::UpdateAlgorithmType updatetype;
	unsigned int verbose;
	std::vector<double> wolfe_constants;

	// secondary options: set by other options
	enum DualityContainerType dualitytype;
	enum MinimizerFactory::InstanceType type;
};



#endif /* COMMANDLINEOPTIONS_HPP_ */
