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
/*! \file algorithms_for_balanced_parentheses.hpp
\brief algorithms.hpp contains algorithms for balanced parentheses sequences.
\author Simon Gog
*/
#ifndef _ALGORITHMS_FOR_BALANCED_PARENTHESES_HPP
#define _ALGORITHMS_FOR_BALANCED_PARENTHESES_HPP

#include "vectors.hpp" // for bit_vector
#include <stack> // for calculate_pioneers_bitmap method
#include <map>   // for calculate_pioneers_bitmap method
#include<algorithm>

namespace algorithm {

		//! Calculate pioneers as defined in the paper of Geary et al. (CPM 2004)
		/*! \param bp The balanced parentheses sequence for that the pioneers should be calculated.
		*  \param block_size Size of the blocks for which the pioneers should be calulated.
		*  \param pioneer_bitmap Reference to the resulting bit_vector.
		*  \par Time complexity
		*       \f$ \Order{n \log n} \f$, where \f$ n=\f$bp.size()
		*  \par Space complexity
		*       \f$ \Order{2n + min(block\_size, \frac{n}{block\_size} )\cdot \log n } \f$
		*/
		template<class size_type>
		void calculate_pioneers_bitmap(const bit_vector &bp, size_type block_size, bit_vector &pioneer_bitmap) {
			pioneer_bitmap.resize(bp.size());      // resize pioneers bitmap
			util::setZeroBits(pioneer_bitmap);  // initialize bitmap with zeros

			std::stack<size_type> opening_parenthesis;
			size_type blocks = (bp.size() + block_size - 1) / block_size;
			// calculate positions of findclose and findopen pioneers
			for (size_type block_nr = 0; block_nr < blocks; ++block_nr) {
				std::map<size_type, size_type> block_and_position; // for find_open and find_close
				std::map<size_type, size_type> matching_position;  // for find_open and find_close
				for (size_type i = 0, j = block_nr*block_size; i < block_size && j < bp.size(); ++i, ++j) {
					if (bp[j]) {//opening parenthesis
						opening_parenthesis.push(j);
					}
					else {// closing parenthesis
						size_type position = opening_parenthesis.top();
						size_type blockpos = position / block_size;
						opening_parenthesis.pop();
						block_and_position[blockpos] = position;
						matching_position[blockpos] = j; // greatest j is pioneer
					}
				}
				for (typename std::map<size_type, size_type>::const_iterator it = block_and_position.begin(),
					end = block_and_position.end(),
					mit = matching_position.begin(); it != end && it->first != block_nr; ++it, ++mit) {
					// opening and closing pioneers are symmetric
					pioneer_bitmap[it->second] = 1;
					pioneer_bitmap[mit->second] = 1;
				}
			}
			// assert that the sequence is balanced
			assert(opening_parenthesis.empty());
		}


		//! Calculates matches (i.e. find_open for a closing parenthesis and find_close for an opening parenthesis) 
		/*! \param bp A bit_vector representing a balanced parentheses sequence.
		*  \param matches Reference to the result.
		*  \pre bp represents a balanced parentheses sequence.
		*  \par Time complexity
		*       \f$ \Order{n} \f$, where \f$ n=\f$bp.size()
		*  \par Space complexity
		*       \f$ \Order{n + 2n\log n } \f$
		*/
		template<class int_vector>
		void calculate_matches(const bit_vector &bp, int_vector &matches) {
			typedef bit_vector::size_type size_type;
			matches = int_vector(bp.size(), 0, bit_magic::l1BP(bp.size()) + 1);
			std::stack<size_type> opening_parenthesis;
			for (size_type i = 0; i < bp.size(); ++i) {
				if (bp[i]) {// opening parenthesis
					opening_parenthesis.push(i);
				}
				else { // closing parenthesis
					assert(!opening_parenthesis.empty());
					size_type position = opening_parenthesis.top();
					opening_parenthesis.pop();
					matches[i] = position;
					assert(matches[i] == position);
					matches[position] = i;
					assert(matches[position] == i);
				}
			}
			// assert that the sequence is balanced
			assert(opening_parenthesis.empty());
		}

		template<class int_vector>
		void calculate_matches_for_pioneers(const bit_vector &bp, const bit_vector &pioneer_bitmap, int_vector &matches) {
			assert(pioneer_bitmap.size() == bp.size());
			typedef bit_vector::size_type size_type;
			matches = int_vector(pioneer_bitmap.size(), 0, bit_magic::l1BP(bp.size()) + 1);
			std::stack<size_type> opening_parenthesis;
			for (size_type i = 0; i < bp.size(); ++i) {
				if (pioneer_bitmap[i]) {
					if (bp[i]) {// opening parenthesis
						opening_parenthesis.push(i);
					}
					else { // closing parenthesis
						assert(!opening_parenthesis.empty());
						size_type position = opening_parenthesis.top();
						opening_parenthesis.pop();
						matches[i] = position;
						assert(matches[i] == position);
						matches[position] = i;
						assert(matches[position] == i);
					}
				}
			}
			// assert that the sequence is balanced
			assert(opening_parenthesis.empty());
		}


		//! Calculates enclose answers for a balanced parentheses sequence.
		/*! \param bp A bit_vector representing a balanced parentheses sequence.
		*  \param enclose Reference to the result.
		*  \pre bp represents a balanced parentheses sequence.
		*  \par Time complexity
		*       \f$ \Order{n} \f$, where \f$ n=\f$bp.size()
		*  \par Space complexity
		*       \f$ \Order{n + 2n\log n } \f$
		*/
		template<class int_vector>
		void calculate_enclose(const bit_vector &bp, int_vector &enclose) {
			typedef bit_vector::size_type size_type;
			enclose = int_vector(bp.size(), 0, bit_magic::l1BP(bp.size()) + 1);
			std::stack<size_type> opening_parenthesis;
			for (size_type i = 0; i < bp.size(); ++i) {
				if (bp[i]) {// opening parenthesis
					if (!opening_parenthesis.empty()) {
						size_type position = opening_parenthesis.top();
						enclose[i] = position;
						assert(enclose[i] == position);
					}
					else
						enclose[i] = bp.size();
					opening_parenthesis.push(i);
				}
				else { // closing parenthesis
					size_type position = opening_parenthesis.top();
					enclose[i] = position; // find open answer if i is a closing parenthesis
					opening_parenthesis.pop();
				}
			}
			// assert that the sequence is balanced
			assert(opening_parenthesis.empty());
		}

		template<class bp_support>
		bool check_bp_support(const bit_vector &bp, bp_support bp_s) {
			typedef bit_vector::size_type size_type;
			// check access and select
			for (size_type i = 0, excess = 0, ones = 0; i < bp.size(); ++i) {
				if (bp[i]) {
					++excess;
					++ones;
					size_type sel = bp_s.select(ones);
					if (sel != i) {
						std::cerr << "select operation: i=" << ones << " value=" << sel << " expected=" << i << std::endl;
						return false;
					}
				}
				else
					--excess;

				if (bp_s.excess(i) != excess) {
					std::cerr << "excess operation: i=" << i << " value=" << bp_s.excess(i) << " expected=" << excess << std::endl;
					return false;
				}
			}
			// check find_open and find_close
			std::stack<size_type> opening_parenthesis;
			for (size_type i = 0; i < bp.size(); ++i) {
				//		std::cerr<<bp[i]<<" i="<<i<<std::endl;
				if (bp[i]) {// opening parenthesis
					opening_parenthesis.push(i);
				}
				else { // closing parenthesis
					assert(!opening_parenthesis.empty());
					size_type position = opening_parenthesis.top();
					size_type fc = bp_s.find_close(position);
					if (fc != i) {
						std::cerr << "find_close operation: i=" << position << " value=" << fc << " expected=" << i << std::endl;
						return false;
					}
					size_type fo = bp_s.find_open(i);
					if (fo != position) {
						std::cerr << "find_open operation: i=" << i << " value=" << fo << " expected=" << position << std::endl;
						return false;
					}
					opening_parenthesis.pop();
				}
			}
			if (!opening_parenthesis.empty()) {
				std::cerr << "balanced parenthese sequence is NOT balanced!" << std::endl;
				return false;
			}
			// check enclose
			for (size_type i = 0; i < bp.size(); ++i) {
				if (bp[i]) {// opening parenthesis
					size_type ec = bp_s.enclose(i);
					size_type position = bp.size();
					if (!opening_parenthesis.empty()) {
						position = opening_parenthesis.top();
					}
					opening_parenthesis.push(i);
					if (ec != position) {
						std::cerr << "encolse operation i=" << i << " value=" << ec << " expected=" << position << std::endl;
					}
				}
				else { // closing parenthesis
					opening_parenthesis.pop();
				}
			}
			return true;
		}

		//! Find the near closing parenthesis if it exists.
		/*!
		* \param bp bit_vector containing the representation of the balanced parentheses sequence.
		* \param i  Position of the opening parenthesis we for which search the corresponding closing parenthesis.
		* \param block_size Number of entries to search for the corresponding closing parenthesis.
		* \return i if there is no near find_close answer, otherwise the position of the near closing parenthesis.
		* \pre We assert that \f$ bp[i]=1 \f$ holds, i.e. there is an opening parenthesis at position i.
		*/
		// TODO: implement a fast version using lookup-tables of size 8
		inline bit_vector::size_type near_find_close_naive(const bit_vector &bp, bit_vector::size_type i, const bit_vector::size_type block_size) {
			typedef bit_vector::size_type size_type;
			size_type opening_parentheses = 1;
			for (size_type j = i + 1; j < (i / block_size + 1)*block_size; ++j) {
				if (bp[j])
					++opening_parentheses;
				else {
					--opening_parentheses;
					if (opening_parentheses == 0) {
						return j;
					}
				}
			}
			return i;
		}

		const uint32_t near_find_close_8bits_lookup[256] = {
			1985229328,2574668850,2574668848,2576971348,2574668816,2576971348,2576971344,2576980342,
			2574668304,2576971346,2576971344,2576980342,2576971280,2576980342,2576980336,2576980377,
			2574660112,2576971314,2576971312,2576980342,2576971280,2576980342,2576980336,2576980377,
			2576970256,2576980338,2576980336,2576980377,2576980240,2576980377,2576980368,2576980377,
			2574529040,2576970802,2576970800,2576980340,2576970768,2576980340,2576980336,2576980377,
			2576970256,2576980338,2576980336,2576980377,2576980240,2576980377,2576980368,2576980377,
			2576953872,2576980274,2576980272,2576980377,2576980240,2576980377,2576980368,2576980377,
			2576978448,2576980370,2576980368,2576980377,2576980240,2576980377,2576980368,2576980377,
			2572431888,2576962610,2576962608,2576980308,2576962576,2576980308,2576980304,2576980377,
			2576962064,2576980306,2576980304,2576980377,2576980240,2576980377,2576980368,2576980377,
			2576953872,2576980274,2576980272,2576980377,2576980240,2576980377,2576980368,2576980377,
			2576978448,2576980370,2576980368,2576980377,2576980240,2576980377,2576980368,2576980377,
			2576626192,2576978994,2576978992,2576980372,2576978960,2576980372,2576980368,2576980377,
			2576978448,2576980370,2576980368,2576980377,2576980240,2576980377,2576980368,2576980377,
			2576953872,2576980274,2576980272,2576980377,2576980240,2576980377,2576980368,2576980377,
			2576978448,2576980370,2576980368,2576980377,2576980240,2576980377,2576980368,2576980377,
			2522100240,2576766002,2576766000,2576979540,2576765968,2576979540,2576979536,2576980374,
			2576765456,2576979538,2576979536,2576980374,2576979472,2576980374,2576980368,2576980377,
			2576757264,2576979506,2576979504,2576980374,2576979472,2576980374,2576980368,2576980377,
			2576978448,2576980370,2576980368,2576980377,2576980240,2576980377,2576980368,2576980377,
			2576626192,2576978994,2576978992,2576980372,2576978960,2576980372,2576980368,2576980377,
			2576978448,2576980370,2576980368,2576980377,2576980240,2576980377,2576980368,2576980377,
			2576953872,2576980274,2576980272,2576980377,2576980240,2576980377,2576980368,2576980377,
			2576978448,2576980370,2576980368,2576980377,2576980240,2576980377,2576980368,2576980377,
			2572431888,2576962610,2576962608,2576980308,2576962576,2576980308,2576980304,2576980377,
			2576962064,2576980306,2576980304,2576980377,2576980240,2576980377,2576980368,2576980377,
			2576953872,2576980274,2576980272,2576980377,2576980240,2576980377,2576980368,2576980377,
			2576978448,2576980370,2576980368,2576980377,2576980240,2576980377,2576980368,2576980377,
			2576626192,2576978994,2576978992,2576980372,2576978960,2576980372,2576980368,2576980377,
			2576978448,2576980370,2576980368,2576980377,2576980240,2576980377,2576980368,2576980377,
			2576953872,2576980274,2576980272,2576980377,2576980240,2576980377,2576980368,2576980377,
			2576978448,2576980370,2576980368,2576980377,2576980240,2576980377,2576980368,2576980377
		};

		const uint32_t near_find_open_8bits_lookup[256] = {
			2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,
			2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,
			2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,
			2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,2576980369,2576980225,
			2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,
			2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,2576980369,2576980225,
			2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,2576980369,2576980225,
			2576980371,2576980371,2576980371,2576980227,2576980259,2576980259,2576978211,2576941347,
			2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,
			2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,2576980369,2576980225,
			2576980377,2576980377,2576980377,2576980377,2576980377,2576980377,2576980369,2576980225,
			2576980371,2576980371,2576980371,2576980227,2576980259,2576980259,2576978211,2576941347,
			2576980373,2576980373,2576980373,2576980373,2576980373,2576980373,2576980373,2576980229,
			2576980373,2576980373,2576980373,2576980229,2576980261,2576980261,2576978213,2576941349,
			2576980293,2576980293,2576980293,2576980293,2576980293,2576980293,2576978245,2576941381,
			2576978757,2576978757,2576978757,2576941893,2576950085,2576950085,2576425797,2566988613,
			2576980375,2576980375,2576980375,2576980375,2576980375,2576980375,2576980375,2576980375,
			2576980375,2576980375,2576980375,2576980375,2576980375,2576980375,2576980375,2576980231,
			2576980375,2576980375,2576980375,2576980375,2576980375,2576980375,2576980375,2576980231,
			2576980375,2576980375,2576980375,2576980231,2576980263,2576980263,2576978215,2576941351,
			2576980375,2576980375,2576980375,2576980375,2576980375,2576980375,2576980375,2576980231,
			2576980375,2576980375,2576980375,2576980231,2576980263,2576980263,2576978215,2576941351,
			2576980295,2576980295,2576980295,2576980295,2576980295,2576980295,2576978247,2576941383,
			2576978759,2576978759,2576978759,2576941895,2576950087,2576950087,2576425799,2566988615,
			2576980327,2576980327,2576980327,2576980327,2576980327,2576980327,2576980327,2576980327,
			2576980327,2576980327,2576980327,2576980327,2576980327,2576980327,2576978279,2576941415,
			2576980327,2576980327,2576980327,2576980327,2576980327,2576980327,2576978279,2576941415,
			2576978791,2576978791,2576978791,2576941927,2576950119,2576950119,2576425831,2566988647,
			2576979303,2576979303,2576979303,2576979303,2576979303,2576979303,2576979303,2576942439,
			2576979303,2576979303,2576979303,2576942439,2576950631,2576950631,2576426343,2566989159,
			2576958823,2576958823,2576958823,2576958823,2576958823,2576958823,2576434535,2566997351,
			2576565607,2576565607,2576565607,2567128423,2569225575,2569225575,2435007847,19088743
		};

		const int8_t excess_8bits_lookup[256] = {
			-8,-6,-6,-4,-6,-4,-4,-2,-6,-4,-4,-2,-4,-2,-2,0,
			-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
			-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
			-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
			-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
			-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
			-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
			-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
			-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
			-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
			-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
			-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
			-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
			-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
			-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
			0,2,2,4,2,4,4,6,2,4,4,6,4,6,6,8
		};

		inline bit_vector::size_type near_find_close(const bit_vector &bp, const bit_vector::size_type i, const bit_vector::size_type block_size) {
			typedef bit_vector::size_type size_type;
			typedef bit_vector::difference_type difference_type;
			difference_type excess = 1;

			const size_type end = ((i + 1) / block_size + 1)*block_size;
			const size_type l = (((i + 1) + 7) / 8) * 8;
			const size_type r = (end / 8) * 8;
			for (size_type j = i + 1; j < std::min(end, l); ++j) {
				if (bp[j])
					++excess;
				else {
					--excess;
					if (excess == 0) {
						return j;
					}
				}
			}
			const uint64_t *b = bp.data();
			for (size_type j = l; j<r; j += 8) {
				if (excess <= 8) {
					assert(excess>0);
					uint32_t x = near_find_close_8bits_lookup[((*(b + (j >> 6))) >> (j & 0x3F)) & 0xFF];
					uint8_t p = (x >> ((excess - 1) << 2)) & 0xF;
					if (p < 9) {
						return j + p;
					}
				}
				excess += excess_8bits_lookup[((*(b + (j >> 6))) >> (j & 0x3F)) & 0xFF];
			}
			for (size_type j = std::max(l, r); j < end; ++j) {
				if (bp[j])
					++excess;
				else {
					--excess;
					if (excess == 0) {
						return j;
					}
				}
			}
			return i;
		}

		// TODO: umbenennen der method in near_find_first_closing_with_excess_diff
		inline bit_vector::size_type near_find_closing(const bit_vector &bp, bit_vector::size_type i, bit_vector::size_type closings, const bit_vector::size_type block_size) {
			typedef bit_vector::size_type size_type;
			typedef bit_vector::difference_type difference_type;
			difference_type excess = 0;
			difference_type succ_excess = -(difference_type)closings;

			const size_type end = (i / block_size + 1)*block_size;
			const size_type l = (((i)+7) / 8) * 8;
			const size_type r = (end / 8) * 8;
			for (size_type j = i; j < std::min(end, l); ++j) {
				if (bp[j])
					++excess;
				else {
					--excess;
					if (excess == succ_excess) {
						return j;
					}
				}
			}
			const uint64_t *b = bp.data();
			for (size_type j = l; j<r; j += 8) {
				if (excess - succ_excess <= 8) {
					//			assert(excess>0);
					uint32_t x = near_find_close_8bits_lookup[((*(b + (j >> 6))) >> (j & 0x3F)) & 0xFF];
					uint8_t p = (x >> (((excess - succ_excess) - 1) << 2)) & 0xF;
					if (p < 9) {
						return j + p;
					}
				}
				excess += excess_8bits_lookup[((*(b + (j >> 6))) >> (j & 0x3F)) & 0xFF];
			}
			for (size_type j = std::max(l, r); j < end; ++j) {
				if (bp[j])
					++excess;
				else {
					--excess;
					if (excess == succ_excess) {
						return j;
					}
				}
			}
			return i - 1;
		}


		//! Find the near opening parenthesis if it exists.
		/*!
		* \param bp bit_vector containing the representation of the balanced parentheses sequence.
		* \param i  Position of the closing parenthesis for which we search the corresponding opening parenthesis.
		* \param block_size Number of entries to search for the corresponding opening parenthesis.
		* \return i if there is no near opening parenthesis, otherwise the position of the near opening parenthesis.
		* \pre We assert that \f$ bp[i]=0 \f$ holds, i.e. there is an closing parenthesis at position i.
		*/
		// TODO: implement a fast version using lookup-tables of size 8
		inline bit_vector::size_type near_find_open_naive(const bit_vector &bp, bit_vector::size_type i, const bit_vector::size_type block_size) {
			typedef bit_vector::size_type size_type;
			size_type closing_parentheses = 1;
			for (size_type j = i; j + block_size - 1 > i && j>0; --j) {
				if (bp[j - 1]) {
					--closing_parentheses;
					if (closing_parentheses == 0) {
						return j - 1;
					}
				}
				else
					++closing_parentheses;
			}
			return i;
		}


		inline bit_vector::size_type near_find_open(const bit_vector &bp, bit_vector::size_type i, const bit_vector::size_type block_size) {
			typedef bit_vector::size_type size_type;
			typedef bit_vector::difference_type difference_type;
			difference_type excess = -1;
			const difference_type begin = ((difference_type)(i - 1) / block_size)*block_size;
			const difference_type r = ((difference_type)(i - 1) / 8) * 8;
			const difference_type l = ((difference_type)((begin + 7) / 8)) * 8;
			for (difference_type j = i - 1; j >= std::max(r, begin); --j) {
				if (bp[j]) {
					if (++excess == 0) {
						return j;
					}
				}
				else
					--excess;
			}
			const uint64_t *b = bp.data();
			for (difference_type j = r - 8; j >= l; j -= 8) {
				if (excess >= -8) {
					assert(excess<0);
					uint32_t x = near_find_open_8bits_lookup[((*(b + (j >> 6))) >> (j & 0x3F)) & 0xFF];
					uint8_t p = (x >> ((-excess - 1) << 2)) & 0xF;
					if (p < 9) {
						return j + p;
					}
				}
				excess += excess_8bits_lookup[((*(b + (j >> 6))) >> (j & 0x3F)) & 0xFF];
			}
			for (difference_type j = std::min(l, r) - 1; j >= begin; --j) {
				if (bp[j]) {
					if (++excess == 0) {
						return j;
					}
				}
				else
					--excess;
			}
			return i;
		}


		inline bit_vector::size_type near_find_opening(const bit_vector &bp, bit_vector::size_type i, const bit_vector::size_type openings, const bit_vector::size_type block_size) {
			typedef bit_vector::size_type size_type;
			typedef bit_vector::difference_type difference_type;
			difference_type excess = 0;
			difference_type succ_excess = openings;

			const difference_type begin = ((difference_type)(i) / block_size)*block_size;
			const difference_type r = ((difference_type)(i) / 8) * 8;
			const difference_type l = ((difference_type)((begin + 7) / 8)) * 8;
			//std::cout<<"i="<<i<<" begin="<<begin<<" r="<<r<<" l="<<l<<std::endl;	
			for (difference_type j = i; j >= std::max(r, begin); --j) {
				if (bp[j]) {
					if (++excess == succ_excess) {
						return j;
					}
				}
				else
					--excess;
			}
			const uint64_t *b = bp.data();
			for (difference_type j = r - 8; j >= l; j -= 8) {
				if (succ_excess - excess <= 8) {
					//			assert(excess<0);
					assert(succ_excess - excess>0);
					uint32_t x = near_find_open_8bits_lookup[((*(b + (j >> 6))) >> (j & 0x3F)) & 0xFF];
					uint8_t p = (x >> ((succ_excess - excess - 1) << 2)) & 0xF;
					if (p < 9) {
						return j + p;
					}
				}
				excess += excess_8bits_lookup[((*(b + (j >> 6))) >> (j & 0x3F)) & 0xFF];
			}
			for (difference_type j = std::min(l, r) - 1; j >= begin; --j) {
				if (bp[j]) {
					if (++excess == succ_excess) {
						return j;
					}
				}
				else
					--excess;
			}
			return i + 1;
		}




		//! Find the opening parenthesis of the enclosing pair if this parenthesis is near.
		/*!
		* \param bp bit_vector containing the representation of the balanced parentheses sequence.
		* \param i Position of the opening parenthesis for which we search the position of the opening parenthesis of the enclosing parentheses pair.
		* \param block_size Number of entries to search for the corresponding opening parenthesis of the enclosing parentheses pair.
		* \return If no near enclose exists return i, otherwise the position of the opening parenthesis of the enclosing pair.
		* \pre We assert that \f$ bp[i]=1 \f$
		*/
		// TODO: implement a fast version using lookup-tables of size 8
		inline bit_vector::size_type near_enclose(const bit_vector &bp, bit_vector::size_type i, const bit_vector::size_type block_size) {
			typedef bit_vector::size_type size_type;
			size_type opening_parentheses = 1;
			for (size_type j = i; j + block_size - 1 > i && j>0; --j) {
				if (bp[j - 1]) {
					++opening_parentheses;
					if (opening_parentheses == 2) {
						return j - 1;
					}
				}
				else
					--opening_parentheses;
			}
			return i;
		}


		const uint16_t min_excess_info_8bits_lookup[256] = {
			17,4105,4360,8201,4615,8713,8456,12297,4870,8968,8968,12297,8711,12809,12552,16393,
			5125,9223,9223,13321,9223,13321,12552,16393,8966,13064,13064,16393,12807,16905,16648,20489,
			5380,9478,9478,13576,9478,13576,13576,16393,9478,13576,13576,16393,12807,16905,16648,20489,
			9221,13319,13319,17417,13319,17417,16648,20489,13062,17160,17160,20489,16903,21001,20744,24585,
			5635,9733,9733,13831,9733,13831,13831,17929,9733,13831,13831,17929,13831,17929,16648,20489,
			9733,13831,13831,17929,13831,17929,16648,20489,13062,17160,17160,20489,16903,21001,20744,24585,
			9476,13574,13574,17672,13574,17672,17672,20489,13574,17672,17672,20489,16903,21001,20744,24585,
			13317,17415,17415,21513,17415,21513,20744,24585,17158,21256,21256,24585,20999,25097,24840,28681,
			5890,9988,9988,14086,9988,14086,14086,18184,9988,14086,14086,18184,14086,18184,18184,20489,
			9988,14086,14086,18184,14086,18184,18184,20489,14086,18184,18184,20489,16903,21001,20744,24585,
			9988,14086,14086,18184,14086,18184,18184,20489,14086,18184,18184,20489,16903,21001,20744,24585,
			13317,17415,17415,21513,17415,21513,20744,24585,17158,21256,21256,24585,20999,25097,24840,28681,
			9731,13829,13829,17927,13829,17927,17927,22025,13829,17927,17927,22025,17927,22025,20744,24585,
			13829,17927,17927,22025,17927,22025,20744,24585,17158,21256,21256,24585,20999,25097,24840,28681,
			13572,17670,17670,21768,17670,21768,21768,24585,17670,21768,21768,24585,20999,25097,24840,28681,
			17413,21511,21511,25609,21511,25609,24840,28681,21254,25352,25352,28681,25095,29193,28936,32777
		};

		inline bit_vector::size_type near_min_excess_position(const bit_vector &bp, const bit_vector::size_type begin, const bit_vector::size_type end) {
			typedef bit_vector::size_type size_type;
			typedef bit_vector::difference_type difference_type;
			difference_type min_excess = end - begin + 1, ex = 0;
			size_type result = end;

			const size_type l = ((begin + 7) / 8) * 8;
			const size_type r = (end / 8) * 8;

			for (size_type k = begin; k < std::min(end, l); ++k) {
				if (bp[k]) {
					++ex;
					if (ex <= min_excess) {
						result = k;
						min_excess = ex;
					}
				}
				else {
					--ex;
				}
			}
			const uint64_t *b = bp.data();// + (l>>6);
			for (size_type k = l; k < r; k += 8) {
				uint16_t x = min_excess_info_8bits_lookup[((*(b + (k >> 6))) >> (k & 0x3F)) & 0xFF];
				int8_t ones = (x >> 12);
				if (ones) {
					int8_t min_ex = (x & 0xFF) - 8;
					if (ex + min_ex <= min_excess) {
						result = k + ((x >> 8) & 0xF);
						min_excess = ex + min_ex;
					}
				}
				ex += ((ones << 1) - 8);
			}
			for (size_type k = std::max(r, l); k < end; ++k) {
				if (bp[k]) {
					++ex;
					if (ex <= min_excess) {
						result = k;
						min_excess = ex;
					}
				}
				else {
					--ex;
				}
			}
			if (min_excess <= ex)
				return result;
			return end;
		}

		inline bit_vector::size_type near_min_excess_position_naive(const bit_vector &bp, const bit_vector::size_type begin, const bit_vector::size_type end) {
			typedef bit_vector::size_type size_type;
			typedef bit_vector::difference_type difference_type;
			difference_type min_excess = end - begin + 1, ex = 0;
			size_type result = end;
			for (size_type k = begin; k<end; ++k) {
				if (bp[k]) {
					++ex;
					if (ex <= min_excess) {
						result = k;
						min_excess = ex;
					}
				}
				else {
					--ex;
				}
			}
			if (min_excess <= ex)
				return result;
			return end;
		}


}// end namespace algorithm

#endif


