#ifndef RMQ_FBST_H
#define RMQ_FBST_H

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
	* Implementation of the range minimum queries based on the binary search tree.
	*
	*  - O(logn) queries
	*  - O(n) preprocessing time
	*  - O(n) words of space
	*
	* In terms of query time, this solution is not as fast as many other solutions (e.g. Fischer11 or Bender05), since
	* it has logarithmic query time and is based on traversing the binary tree at most two times. On the positive side,
	* it requires linear space and has good maintainability due the construction simplicity. The construction seems to be
	* the fastest one among the most popular methods. Despite of slower query times, this approach is the most economical one.
	*
	*/
	template<typename T>
	class rmq_fbst {
	public:

		rmq_fbst(std::vector<T>& v) : n(v.size()), a(v) {}

		//----------------------------------------------------------------------------------------------------------------

		int operator[] (int pos) { return a[pos]; }

		//----------------------------------------------------------------------------------------------------------------

		void init_sequential()
		{
			q.push_back(std::vector<int>());

			for (int i = 1; i < n; i += 2)
			{
				if (a[i - 1] < a[i]) { q[0].push_back(i - 1); }
				else { q[0].push_back(i); }
			}

			for (int i = 0; i < std::ceil(std::log2(n)) - 1; ++i)
			{
				q.push_back(std::vector<int>());

				std::vector<int>& last = q[q.size() - 2];
				std::vector<int>& curr = q[q.size() - 1];

				for (int i = 1; i < last.size(); i += 2)
				{
					if (a[last[i - 1]] < a[last[i]]) { curr.push_back(last[i - 1]); }
					else { curr.push_back(last[i]); }
				}

				if (last.size() & 2 != 0) { curr.push_back(last[last.size() - 1]); }
			}
		}

		//----------------------------------------------------------------------------------------------------------------

		int get_minimum(int i, int j)
		{
			int min1 = std::numeric_limits<int>::max();
			int min2 = std::numeric_limits<int>::max();
			int min3 = std::numeric_limits<int>::max();

			int delta = j - i + 1;
			int level = std::floor(std::log2(delta));
			int range = util::ipow(2, level);

			int b1 = i / range;
			int b2 = j / range;

			if (i > b1*range) { min1 = get_minimum(i, b1*range + range - 1); b1++; }
			if (j < b2*range + range - 1) { min3 = get_minimum(b2*range, j); b2--; }

			if (b1 <= b2)
			{
				int first = (level == 0) ? b1 : q[level - 1][b1];
				int second = (level == 0) ? b2 : q[level - 1][b2];
				if (a[first] < a[second]) { min2 = first; }
				else { min2 = second; }
			}
			else { if (a[i] < a[j]) { min2 = i; } else { min2 = j; } }

			return get_min(min1, min2, min3);
		}

		//----------------------------------------------------------------------------------------------------------------

		int size() { return n; }
		int size_in_bits() { int count = 0; for (int i = 0; i < q.size(); ++i) count += q[i].size(); return count * sizeof(int) * 8; }

		//----------------------------------------------------------------------------------------------------------------

	private:

		int get_min(int x, int y, int z) { return get_min(x, get_min(y, z)); }

		int get_min(int x, int y)
		{
			if (x == std::numeric_limits<int>::max()) return y;
			else if (y == std::numeric_limits<int>::max()) return x;
			else if (a[x] < a[y]) { return x; }
			else { return y; }
		}

		//----------------------------------------------------------------------------------------------------------------

		std::vector<T>& a; // input
		int n; // input size
		std::vector<std::vector<int>> q; // RMQ structure
	};

};

#endif // RMQ_FBST_H
