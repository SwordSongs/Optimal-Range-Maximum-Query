/* sdsl - succinct data structures library
Copyright (C) 2009 Simon Gog

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file algorithms_for_compressed_suffix_trees.hpp
\brief algorithms_for_compressed_suffix_trees.hpp contains algorithms for compressed suffix trees.
\author Simon Gog
*/
#ifndef _ALGORITHMS_FOR_COMPRESSED_SUFFIX_TREES_HPP
#define _ALGORITHMS_FOR_COMPRESSED_SUFFIX_TREES_HPP

#include "vectors.hpp" // for bit_vector
#include <stack> // for calculate_supercartesian_tree_bp

namespace algorithm {

		//! Calculate the balanced parentheses of the Super-Cartesian tree, described in Ohlebusch and Gog (SPIRE 2009). 
		/*! \param vec Random access container for which the Super-Cartesian tree representation should be calculated.
		*             The value_type of vec should be an unsigned integer type.
		*  \param bp Reference to the balanced parentheses sequence which represents the Super-Cartesian tree.
		*  \par Time complexity
		*       \f$ \Order{2n} \f$, where \f$ n=\f$vec.size()
		*  \par Space complexity
		*       \f$ \Order{n \cdot \log n } \f$ bits, can be improved to \f$\Order{n}\f$ bits, by storing the differences on a succinct stack.
		*/
		template<class RandomAccessContainer>
		void construct_supercartesian_tree_bp(const RandomAccessContainer &vec, bit_vector &bp) {
			typedef typename RandomAccessContainer::size_type size_type;
			bp.resize(2 * vec.size());    // resize bit vector for balanaced parantheses to 2 n bits
			std::stack<typename RandomAccessContainer::value_type> vec_stack;
			vec_stack.push(0);
			bp[0] = 1; // writing an opening parenthesis
			size_type k = 1;
			for (size_type i = 1, l; i < vec.size(); ++i) {
				l = vec[i];
				while (vec_stack.size() > 0 and l < vec_stack.top()) {
					vec_stack.pop();
					bp[k++] = 0; // writing a closing parenthesis
				}
				vec_stack.push(l);
				bp[k++] = 1; // writing an opening  parenthesis
			}
			while (vec_stack.size() > 0) {
				vec_stack.pop();
				bp[k++] = 0; // writing a closing parenthesis
			}
			assert(k == 2 * vec.size());
		}

}// end namespace algorithm

#endif

