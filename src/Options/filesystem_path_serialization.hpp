/*
 * filesystem_path_serialization.hpp
 *
 *  Created on: Jun 25, 2015
 *      Author: heber
 */

#ifndef OPTIONS_FILESYSTEM_PATH_SERIALIZATION_HPP_
#define OPTIONS_FILESYSTEM_PATH_SERIALIZATION_HPP_

#include "BassoConfig.h"

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <boost/config.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/level.hpp>
#include <boost/serialization/nvp.hpp>

//BOOST_CLASS_IMPLEMENTATION(boost::filesystem::path, boost::serialization::object_serializable)
//#ifndef BOOST_NO_STD_WSTRING
//BOOST_CLASS_IMPLEMENTATION(boost::filesystem::wpath, boost::serialization::object_serializable)
//#endif

namespace boost { namespace serialization {

template<class Archive>
void serialize(
		Archive& ar,
		boost::filesystem::path& p,
		const unsigned int version)
{
     std::string s;
     if(Archive::is_saving::value)
         s = p.string();
     ar & boost::serialization::make_nvp("string", s);
     if(Archive::is_loading::value)
         p = s;
}

} /* namespace serialization */

} /* namespace boost */

#endif /* OPTIONS_FILESYSTEM_PATH_SERIALIZATION_HPP_ */
