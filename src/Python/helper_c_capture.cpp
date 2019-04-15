/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2019 Frederik Heber. All rights reserved.
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

/* helper_c_capture.cpp
 *
 * This is inspired by the answer on stackoverflow
 * https://stackoverflow.com/questions/35745541/how-to-get-printed-output-from-ctypes-c-functions-into-jupyter-ipython-notebook
 * where certain libc functions
 * are exposed using a python module (with ctypes). This is however very hacky
 * and as I am using a C++ module anyway, I rather expose these functions cleanly
 * through an extra library.
 *
 * Note:
 *   It was not straight-forward to expose _IO_FILE to boost::python. Hence,
 *   I simply crippled the exposed function to have not arguments and work on
 *   stdout directly.
 *
 * See \b create_c_contextmanager.py
 *
 *  Created on: Apr 15, 2019
 *      Author: heber
 */

#include <boost/python.hpp>

#include <cstdio>

using namespace boost::python;

void fflush_stdout() { fflush(stdout); }
void disable_stdout_buffering() { setbuf(stdout, NULL); }
void fflush_stderr() { fflush(stderr); }
void disable_stderr_buffering() { setbuf(stderr, NULL); }

BOOST_PYTHON_MODULE(basso_capture)
{
	def("fflush_stderr", fflush_stderr, "flushing stderr buffer");
	def("fflush_stdout", fflush_stdout, "flushing stdout buffer");
	def("disable_stderr_buffering", disable_stderr_buffering,
			"disabling the error buffer");
	def("disable_stdout_buffering", disable_stdout_buffering,
			"disabling the output buffer");
}
