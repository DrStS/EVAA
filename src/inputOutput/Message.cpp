/*
 * Copyright &copy; 2019, Dr. Stefan Sicklinger, Munich \\n
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

#include "Message.h"
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/phoenix/core/expression.hpp>
#include <boost/phoenix/bind/bind_function.hpp>
#include <boost/log/utility/setup/file.hpp>

namespace EVAA {

std::string_view formatMessage(
    boost::log::value_ref<
        std::string, boost::log::expressions::tag::smessage> const &t_message) {
  // Check to see if the attribute value has been found
  if (t_message) {
    std::string_view msg = t_message.get();
    if (!msg.empty() && msg.back() == '\n') {
      msg = std::string_view(msg.data(), msg.size() - 1);
    }
    return msg;
  }

  return std::string_view();
}

void Message::initLogging() {

  boost::log::add_file_log(boost::log::keywords::file_name = "EVAA.log",
                           boost::log::keywords::target_file_name = "EVAA.log");
  boost::log::add_console_log(
      std::clog,
      boost::log::keywords::format =
          (boost::log::expressions::stream
           << boost::log::expressions::format_date_time<
                  boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S")
           << " ["
           << boost::log::expressions::attr<
                  boost::log::trivial::severity_level>("Severity")
           << "]: "
           << boost::phoenix::bind(&formatMessage,
                                   boost::log::expressions::smessage // NOLINT
                                       .or_none())));
  boost::log::add_common_attributes();
  // Set filter based on XML input
  boost::log::core::get()->set_filter(boost::log::trivial::severity >=
                                      boost::log::trivial::debug);
}
// clang-format off
	void Message::writeASCIIArt() {
		LOG_INFO << "#            _____                    _____                    _____                    _____                              " << std::endl;
		LOG_INFO << R"(#           /\    \                  /\    \                  /\    \                  /\    \                     )" << std::endl;
		LOG_INFO << R"(#          /::\    \                /::\____\                /::\    \                /::\    \                    )" << std::endl;
		LOG_INFO << R"(#         /::::\    \              /:::/    /               /::::\    \              /::::\    \                     )" << std::endl;
		LOG_INFO << R"(#        /::::::\    \            /:::/    /               /::::::\    \            /::::::\    \                    )" << std::endl;
		LOG_INFO << R"(#       /:::/\:::\    \          /:::/    /               /:::/\:::\    \          /:::/\:::\    \                )" << std::endl;
		LOG_INFO << R"(#      /:::/__\:::\    \        /:::/____/               /:::/__\:::\    \        /:::/__\:::\    \               )" << std::endl;
		LOG_INFO << R"(#     /::::\   \:::\    \       |::|    |               /::::\   \:::\    \      /::::\   \:::\    \           )" << std::endl;
		LOG_INFO << R"(#    /::::::\   \:::\    \      |::|    |     _____    /::::::\   \:::\    \    /::::::\   \:::\    \          )" << std::endl;
		LOG_INFO << R"(#   /:::/\:::\   \:::\    \     |::|    |    /\    \  /:::/\:::\   \:::\    \  /:::/\:::\   \:::\    \    )" << std::endl;
		LOG_INFO << R"(#  /:::/__\:::\   \:::\____\    |::|    |   /::\____\/:::/  \:::\   \:::\____\/:::/  \:::\   \:::\____\   )" << std::endl;
		LOG_INFO << R"(#  \:::\   \:::\   \::/    /    |::|    |  /:::/    /\::/    \:::\  /:::/    /\::/    \:::\  /:::/    /         )" << std::endl;
		LOG_INFO << R"(#   \:::\   \:::\   \/____/     |::|    | /:::/    /  \/____/ \:::\/:::/    /  \/____/ \:::\/:::/    /          )" << std::endl;
		LOG_INFO << R"(#    \:::\   \:::\    \         |::|____|/:::/    /            \::::::/    /            \::::::/    /               )" << std::endl;
		LOG_INFO << R"(#     \:::\   \:::\____\        |:::::::::::/    /              \::::/    /              \::::/    /                )" << std::endl;
		LOG_INFO << R"(#      \:::\   \::/    /        \::::::::::/____/               /:::/    /               /:::/    /                    )" << std::endl;
		LOG_INFO << R"(#       \:::\   \/____/          ~~~~~~~~~~                    /:::/    /               /:::/    /                      )" << std::endl;
		LOG_INFO << R"(#        \:::\    \                                           /:::/    /               /:::/    /                       )" << std::endl;
		LOG_INFO << R"(#         \:::\____\                                         /:::/    /               /:::/    /                        )" << std::endl;
		LOG_INFO << R"(#          \::/    /                                         \::/    /                \::/    /                         )" << std::endl;
		LOG_INFO << R"(#           \/____/                                           \/____/                  \/____/                          )" << std::endl;
		LOG_INFO << "#                                                                                                                          " << std::endl;
	};
// clang-format on
} // namespace EVAA
