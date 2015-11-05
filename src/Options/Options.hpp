/*
 * Options.hpp
 *
 *  Created on: Jul 27, 2015
 *      Author: heber
 */

#ifndef OPTIONS_OPTIONS_HPP_
#define OPTIONS_OPTIONS_HPP_

#include "BassoConfig.h"

#include <iosfwd>
#include <string>

#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>

#include "Options/filesystem_path_serialization.hpp"

/** This class defines the interface for all derived instances that obtain
 * options from the command-line options and store them internally.
 *
 */
class Options
{
public:
	Options();

	virtual ~Options() {}

	/** Initializes the options.
	 *
	 */
	virtual void init();

	/** Performs the actual parsing.
	 *
	 * \param argc argument count
	 * \param **argv array of arguments
	 */
	virtual void parse(int argc, char **argv);

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
	virtual bool checkSensibility() const { return true; }

	/** Sets other values that require checkSensibility().
	 *
	 */
	virtual void setSecondaryValues() {}

	/** Stores all variables as "key = value" pairs per line under the given
	 * stream \a _output.
	 *
	 * @param _output output stream to write to
	 */
	virtual void store(std::ostream &_output) const;


protected:
	/** Helper template to write "key = value" pair to a stream.
	 *
	 * Note that the template class \a T designates the type of the
	 * value.
	 *
	 * \note This function resides here to allow use by derived classes
	 * while not requiring to have a new namespace (i.e. we just use
	 * Options as a namespace).
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
	boost::filesystem::path config_filename;
	//!> desired verbosity level
	unsigned int verbose;

private:
	friend class boost::serialization::access;
	template<class Archive>
	void save(Archive & ar, const unsigned int version) const
	{
		ar & config_filename;
		ar & verbose;
	}

	template<class Archive>
	void load(Archive & ar, const unsigned int version)
	{
		ar & config_filename;
		ar & verbose;

		// set verbosity after loading
		setVerbosity();
	}

	BOOST_SERIALIZATION_SPLIT_MEMBER()
};

#endif /* OPTIONS_OPTIONS_HPP_ */
