/*
 * Copyright &copy; 2020, Dr. Stefan Sicklinger, Munich \n
 *
 *  All rights reserved.
 *
 *  This file is part of EVAA.
 *
 *  EVAA is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EVAA is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EVAA.  If not, see http://www.gnu.org/licenses/.
 */

#pragma once
#include <string>
#include <boost/log/trivial.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>

namespace EVAA {
/**
 * \brief Class wrapper to boost logging
 **/
class Message {
public:
  /**
   * \brief Write EVAA ASCII Art
   * \author Stefan Sicklinger
   **/
  void writeASCIIArt();
  /**
   * \brief init logger
   * \author Stefan Sicklinger
   **/
  void initLogging();

private:
  /// This variable is shared between all objects of type message
};
std::string_view formatMessage(
    boost::log::value_ref<
        std::string, boost::log::expressions::tag::smessage> const &t_message);

#define LOG_DEBUG                                                              \
  BOOST_LOG_STREAM_WITH_PARAMS(                                                \
      ::boost::log::trivial::logger::get(),                                    \
      (::boost::log::keywords::severity = ::boost::log::trivial::debug))

#define LOG_INFO                                                               \
  BOOST_LOG_STREAM_WITH_PARAMS(                                                \
      ::boost::log::trivial::logger::get(),                                    \
      (::boost::log::keywords::severity = ::boost::log::trivial::info))

#define LOG_WARNING                                                            \
  BOOST_LOG_STREAM_WITH_PARAMS(                                                \
      ::boost::log::trivial::logger::get(),                                    \
      (::boost::log::keywords::severity = ::boost::log::trivial::warning))

#define LOG_ERROR                                                              \
  BOOST_LOG_STREAM_WITH_PARAMS(                                                \
      ::boost::log::trivial::logger::get(),                                    \
      (::boost::log::keywords::severity = ::boost::log::trivial::error))

} // namespace EVAA
