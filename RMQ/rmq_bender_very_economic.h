#ifndef RMQ_BENDER_VERY_ECONOMIC_H
#define RMQ_BENDER_VERY_ECONOMIC_H

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
#include "rmq_bender.h"

namespace util {

	/**
	* Implementation of the range minimum queries according to Bender 05 (very economic variant).
	*
	* We achieve linear space usage (very close to FBST) at cost of logarithmic query times due in-block scanning. The hidden constant
	* in the space complexity is relatively small and the query times are very fast due the small values of logn in practice,
	* so that avoiding one ceil and one log2 calculations for each query turns profitable. In practice, this variant seems to be a potential
	* candidate to any other implementation: one of the fastest query times and very low space consumption.
	*
	*  - O(logn) queries
	*  - O(n) preprocessing time
	*  - O(n) words of space
	*
	*/
	template<typename T>
	class rmq_bender00_very_eco {
	public:

		rmq_bender00_very_eco(std::vector<T>& v)
			: a(v), directory_rmq(directory)
		{
			n = v.size();
			block_size = std::ceil(std::log2(n));
			nblocks = util::idiv_ceil(n, block_size);

			b.resize(nblocks, std::vector<T>());
		}

		//----------------------------------------------------------------------------------------------------------------

		int operator[] (int pos) { return b[pos / block_size][pos % block_size]; }

		//----------------------------------------------------------------------------------------------------------------

		void init_sequential()
		{
			int min_idx = 0;
			int min = std::numeric_limits<int>::max();
			int curr_block = 0;

			for (int i = 0; i < n; ++i)
			{
				b[curr_block].push_back(a[i]);

				if (a[i] < min) { min = a[i]; min_idx = i; }

				if ((i + 1) % block_size == 0)
				{
					curr_block++;
					directory_inv.push_back(min_idx);
					directory.push_back(min);
					min = std::numeric_limits<int>::max();
				}
			}

			directory_rmq.init_sequential();
			directory.clear();
		}

		//----------------------------------------------------------------------------------------------------------------

		int get_minimum(int i, int j)
		{
			util::result res1, res2, res3;

			int b1 = i / block_size;
			int b2 = j / block_size;

			if (i > b1*block_size)
			{
				res1.idx = util::get_minimum_index(i % block_size, (j < (b1 + 1)*block_size) ? (j % block_size) : (block_size - 1), b[b1]);
				res1.val = b[b1][res1.idx];
				res1.idx += b1*block_size;
				b1++;
			}

			if (j < b2*block_size + block_size - 1)
			{
				res3.idx = util::get_minimum_index((i >= b2*block_size) ? (i % block_size) : 0, j % block_size, b[b2]);
				res3.val = b[b2][res3.idx];
				res3.idx += b2*block_size;
				b2--;
			}

			if (b1 <= b2)
			{
				res2.idx = directory_inv[directory_rmq.get_minimum(b1, b2)];
				res2.val = b[res2.idx / block_size][res2.idx % block_size];
			}

			return get_min(res1, res2, res3).idx;
		}

		//----------------------------------------------------------------------------------------------------------------

		int size() { return n; }

		int size_in_bits()
		{
			int inverse_size = directory_inv.size() * sizeof(int) * 8;
			int directory_size = directory_rmq.size_in_bits();

			return inverse_size + directory_size;
		}

		//----------------------------------------------------------------------------------------------------------------

	private:

		util::result get_min(util::result a, util::result b)
		{
			if (a.idx == std::numeric_limits<int>::max()) { return b; }
			else if (b.idx == std::numeric_limits<int>::max()) { return a; }
			else if (a < b) { return a; }
			else { return b; }
		}

		util::result get_min(util::result a, util::result b, util::result c) { return get_min(a, get_min(b, c)); }

		//----------------------------------------------------------------------------------------------------------------

		std::vector<T>& a; // input

						   /** parameters */
		int n; // input size
		int nblocks; // n/logn
		int block_size; // logn

						/** RMQ structure of the blocks of size logn*/
		std::vector<std::vector<int>> b;

		/** RMQ structure of the directory with the minimum element of each block */
		std::vector<int> directory;
		std::vector<int> directory_inv; // n/logn additional words
		util::rmq_bender00<T> directory_rmq;
	};

};

#endif // RMQ_BENDER_VERY_ECONOMIC_H