#ifndef RMQ_FISCHER_H
#define RMQ_FISCHER_H

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
	* Implementation of the range minimum queries according to Fischer 11.
	*
	* The major drawback of this technique is the slow construction time in practice. However, it results in
	* the fast queries, since only memory accesses and few comparisons are required. Moreover, the EXPECTED constant
	* involved in the space complexity is not big, because only one extra integer is stored for every block,
	* and not a complete RMQ structure. This solution is recommended in case the expected space consumption is a very
	* important issue. To overcome the precalculation problem, this step can be easily parallelized.
	*
	*  - O(1) queries
	*  - O(sn) preprocessing time, if s = logn is the block size
	*  - O(n) words of space
	*
	*/
	template<typename T>
	class rmq_fischer11 {
	public:

		rmq_fischer11(std::vector<T>& v)
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

			types.resize(nblocks);
			for (int i = 0; i < n; ++i)
			{
				b[curr_block].push_back(a[i]);

				if (a[i] < min) { min = a[i]; min_idx = i; }

				if ((i + 1) % block_size == 0)
				{
					types[curr_block] = get_type(b[curr_block++]);
					directory_inv.push_back(min_idx);
					directory.push_back(min);
					min = std::numeric_limits<int>::max();
				}
			}

			directory_rmq.init_sequential();
			directory.clear();

			// allocate space for precalculations of the in-block minima
			int max_type = util::get_maximum(0, types.size() - 1, types);

			c.resize(block_size, std::vector<std::vector<int>>());

			for (int i = 0; i < block_size; ++i)
			{
				c[i].resize(block_size, std::vector<int>());

				for (int j = 0; j < block_size; ++j)
				{
					c[i][j].resize(max_type + 1, -1);
				}
			}

			// in-block minima precalculations
			for (int idx = 0; idx < types.size(); ++idx)
			{
				for (int i = 0; i < block_size; ++i)
				{
					for (int j = i; j < block_size; ++j)
					{
						c[i][j][types[idx]] = util::get_minimum_index(i, j, b[idx]);
					}
				}
			}
		}

		//----------------------------------------------------------------------------------------------------------------

		int get_minimum(int i, int j)
		{
			util::result res1, res2, res3;

			int b1 = i / block_size;
			int b2 = j / block_size;

			if (i > b1*block_size)
			{
				res1.idx = c[i % block_size][(j < (b1 + 1)*block_size) ? (j % block_size) : (block_size - 1)][types[b1]];
				res1.val = b[b1][res1.idx];
				res1.idx += b1*block_size;
				b1++;
			}

			if (j < b2*block_size + block_size - 1)
			{
				res3.idx = c[(i >= b2*block_size) ? (i % block_size) : 0][j % block_size][types[b2]];
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
			int lookup_size = block_size*block_size*types[types.size() - 1];
			int types_size = types.size() * sizeof(int) * 8;
			int inverse_size = directory_inv.size() * sizeof(int) * 8;
			int directory_size = directory_rmq.size_in_bits();

			return lookup_size + types_size + inverse_size + directory_size;
		}

		//----------------------------------------------------------------------------------------------------------------

	private:

		int get_type(std::vector<int>& block)
		{
			int s = block.size();

			std::vector<int> r(s + 1, std::numeric_limits<int>::min());
			int q = s; int N = 0;

			for (int i = 0; i < s; i++)
			{
				while (r[q + i - s] > block[i])
				{
					N += ballot_number(s - i, q);
					q--;
				}

				r[q + i + 1 - s] = block[i];
			}

			return N;
		}

		int ballot_number(int p, int q)
		{
			if (p == 0 && q == 0) { return 1; }
			else if (p >= 0 && q >= p && q != 0) { return ballot_number(p, q - 1) + ballot_number(p - 1, q); }
			else { return 0; }
		}

		//----------------------------------------------------------------------------------------------------------------

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
		std::vector<std::vector<std::vector<int>>> c;
		std::vector<int> types;
		std::vector<std::vector<int>> b;

		/** RMQ structure of the directory with the minimum element of each block */
		std::vector<int> directory;
		std::vector<int> directory_inv; // n/logn additional words
		util::rmq_bender00<T> directory_rmq;
	};

};

#endif // RMQ_FISCHER_H