#ifndef RMQ_BENDER_H
#define RMQ_BENDER_H

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

#include <vector>
#include <cmath>
#include "util.h"

namespace util {

	/**
	* Implementation of the range minimum queries according to Bender 05.
	*
	* The most relevant advantage of this technique is the constant query time, whereas it is easy to give maintenence due the simplicity
	* of the involved methods. In comparison to the economic variant of the Bender's technique, the queries are answered much faster in practice
	* due the simple calculations. To perform a query, we need two math calculations, namely a floor and a base two logarithm, two memory access
	* for fetching precalculated queries on ranges and one integer comparison. These operations are very cheap in any modern computer.
	* On the negative side, it is more expensive in terms of space usage (extra logn factor or 3-4x more space for n=10,000 against the
	* economic variant). The construction method is via dynamic programming by solving the recurrence q[i][l] = argmin{q[i][l - 1], q[i + 2^l][l-1]}.
	*
	*  - O(1) queries
	*  - O(nlogn) preprocessing time
	*  - O(nlogn) words of space
	*
	*/
	template<typename T>
	class rmq_bender00 {
	public:

		rmq_bender00(std::vector<T>& v) : n(v.size()), a(v), q(v.size(), std::vector<int>()) {}

		//----------------------------------------------------------------------------------------------------------------

		int operator[] (int pos) { return a[pos]; }

		//----------------------------------------------------------------------------------------------------------------

		void init_sequential()
		{
			n = a.size(); q.resize(n);

			for (int i = 0; i < n; ++i) { q[i].push_back(get_min(i, std::min(i + 1, n - 1))); }

			for (int l = 1; l < std::floor(std::log2(n)); ++l)
			{
				for (int i = 0; i < n; ++i)
				{
					q[i].push_back(get_min(q[i][l - 1], q[std::min(n - 1, i + util::ipow(2, l))][l - 1]));
				}
			}
		}

		//----------------------------------------------------------------------------------------------------------------

		int get_minimum(int i, int j)
		{
			if (i == j) { return i; }
			else {
				int h = std::floor(std::log2(j - i));
				if (h == 0) { return q[i][0]; }
				else { return get_min(q[i][h - 1], q[j - util::ipow(2, h) + 1][h - 1]); }
			}
		}

		//----------------------------------------------------------------------------------------------------------------

		int size() { return n; }

		int size_in_bits() { int count = 0; for (int i = 0; i < q.size(); ++i) { count += q[i].size(); } return count * sizeof(int) * 8; }

		//----------------------------------------------------------------------------------------------------------------

	private:

		int get_min(int i, int j) { return a[i] < a[j] ? i : j; }

		//----------------------------------------------------------------------------------------------------------------

		std::vector<T>& a; // input
		int n; // input size
		std::vector<std::vector<int>> q; // RMQ structure
	};

};

#endif // RMQ_BENDER_H