#ifndef UTIL_H
#define UTIL_H

/*
* Copyright (c) 2016, Fellipe Bernardes Lima <fellipebl@gmail.com>.
* All rights reserved. Redistribution and use in source and binary forms,
* with or without modification, are permitted provided that the following
* conditions are met:
*
* (1) Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
*
* (2) Redistributions in binary form must reproduce the above copyright
* notice, this list of conditions and the following disclaimer in the
* documentation and/or other materials provided with the distribution.
*
* (3) The names of the creators and contributors may not be used to endorse
* or promote products derived from this software without specific prior
* written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*/

#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <ctime>

//----------------------------------------------------------------------------------------------------------------

time_t time_begin;
time_t time_end;
void begin_timing() { time_begin = clock(); }
void end_timing() { time_end = clock(); std::cout << "Elapsed time: " << (time_end - time_begin) / CLOCKS_PER_SEC << " secs!" << std::endl; }

//----------------------------------------------------------------------------------------------------------------

namespace util {

	struct result {
		result()
		{
			idx = std::numeric_limits<int>::max();
			val = std::numeric_limits<int>::max();
		}

		bool operator< (result& rhs) { return (*this).val < rhs.val; }

		int idx;
		int val;
	};

	//----------------------------------------------------------------------------------------------------------------

	// division without overflows
	int idiv_ceil(int dividend, int divisor) { return 1 + (dividend - 1) / divisor; }

	// power of integers via exponentiation by squaring: calculate n^e in O(logn) time
	int ipow(int n, int e) { int r = 1; while (e) { if (e & 1) { r *= n; } n *= n; e >>= 1; } return r; }

	//----------------------------------------------------------------------------------------------------------------

	template<typename S, typename T>
	void get_minima_scanning(int i, int j, S& array, T& res)
	{
		int min = std::numeric_limits<int>::max();
		for (int k = i; k <= j; ++k) { if (array[k] < min) { min = array[k]; } }
		for (int k = i; k <= j; ++k) { if (array[k] == min) { res.push_back(k); } }
	}

	template<typename S>
	int get_minimum(int i, int j, S& array)
	{
		int min = std::numeric_limits<int>::max();
		for (int k = i; k <= j; ++k) { if (array[k] < min) { min = array[k]; } }
		return min;
	}

	template<typename S>
	int get_maximum(int i, int j, S& array)
	{
		int max = std::numeric_limits<int>::min();
		for (int k = i; k <= j; ++k) { if (array[k] > max) { max = array[k]; } }
		return max;
	}

	template<typename S>
	int get_minimum_index(int i, int j, S& array)
	{
		int min = std::numeric_limits<int>::max();
		int min_idx = std::numeric_limits<int>::max();
		for (int k = i; k <= j; ++k) if (array[k] < min) { min_idx = k; min = array[k]; }
		return min_idx;
	}

	//----------------------------------------------------------------------------------------------------------------

	void print_vec(std::vector<int>& in)
	{
		for (int i = 0; i < in.size(); ++i) { std::cout << std::setw(std::log(in.size()) + 1) << i << std::flush; } std::cout << std::endl;
		for (int i = 0; i < in.size(); ++i) { std::cout << std::setw(std::log(in.size()) + 1) << in[i] << std::flush; } std::cout << std::endl;
	}

}

#endif // UTIL_H