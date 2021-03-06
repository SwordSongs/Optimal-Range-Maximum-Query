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
/*! \file elias_delta_coder.hpp
\brief elias_delta_coder.hpp contains the class sdsl::coder::elias_delta
\author Simon Gog
*/
#ifndef _ELIAS_DELTA_CODER_HPP
#define _ELIAS_DELTA_CODER_HPP

#include "vectors.hpp"

namespace coder {

		//! A class to encode and decode between Elias-\f$\delta\f$ and binary code.
		class elias_delta {
		public:
			typedef uint64_t size_type;

			//! Array contains precomputed values for the decoding of the prefix sum of Elias-Delta encoded numbers.
			/*! The 8 most significant bits contain the length of decoded bits.
			*  The following 8 bits contain the number of decoded values.
			*  The last 16 bits contain the sum of the decoded values.
			*/
			static const uint32_t EliasDeltaPrefixSum[1 << 16];

			static const uint16_t EliasDeltaPrefixSum8bit[(1 << 8) * 8];
			static const uint8_t min_codeword_length = 1; // 1 represents 1 and is the code word with minimum length
			static uint8_t encoding_length(uint64_t);
			//! Decode n Elias-delta encoded bits beginning at start_idx in the bitstring "data"
			/* \param data Bitstring
			\param start_idx Starting index of the decoding.
			\param n Number of values to decode from the bitstring.
			\param it Iterator to decode the values.
			*/
			template<bool sumup, bool increment, class Iterator>
			static uint64_t decode(const uint64_t *data, const size_type start_idx, size_type n, Iterator it = (Iterator)NULL);

			//! Decode n EliasDelta encoded integers beginning at start_idx in the bitstring "data"  and return the sum of these values.
			/*! \param data Pointer to the beginning of the EliasDelta encoded bitstring.
			\param start_idx Index of the first bit to endcode the values from.
			\param n Number of values to decode from the bitstring. Attention: There have to be at least n encoded values in the bitstring.
			*/
			static uint64_t decode_prefix_sum(const uint64_t *data, const size_type start_idx, size_type n);
			static uint64_t decode_prefix_sum(const uint64_t *data, const size_type start_idx, const size_type end_idx, size_type n);

			template<class int_vector>
			static bool encode(const int_vector &v, int_vector &z);
			template<class int_vector>
			static bool decode(const int_vector &z, int_vector &v);

			//! Encode one positive integer x to an int_vector at bit position start_idx.
			/* \param x Positive integer to encode.
			\param z Raw data of vector to write the encoded form of x.
			\param start_idx Beginning bit index to write the encoded form ox x in z.
			*/
			static void encode(uint64_t x, uint64_t *&z, uint8_t &offset);

			template<class int_vector>
			static uint64_t* raw_data(int_vector &v) {
				return v.m_data;
			};
		};

		// \sa coder::elias_delta::encoding_length
		inline uint8_t elias_delta::encoding_length(uint64_t w) {
			uint8_t len_1 = bit_magic::l1BP(w);
			return len_1 + (bit_magic::l1BP(len_1 + 1) << 1) + 1;
		}

		template<class int_vector>
		bool elias_delta::encode(const int_vector &v, int_vector &z) {
			z.setIntWidth(v.getIntWidth());
			//	z.m_elements = v.size();
			size_t z_bit_size = 0;
			for (typename int_vector::const_iterator it = v.begin(), end = v.end(); it != end; ++it) {
				z_bit_size += encoding_length(*it);
			}
			z.bit_resize(z_bit_size); // Initial size of z
			if (z_bit_size & 0x3F) { // if z_bit_size % 64 != 0
				*(z.m_data + (z_bit_size >> 6)) = 0; // initialize last word
			}
			z_bit_size = 0;
			uint64_t *z_data = z.m_data;
			uint8_t offset = 0;
			uint64_t w;
			size_t len, len_1_len; // TODO: change to uint8_t and test it
			for (typename int_vector::const_iterator it = v.begin(), end = v.end(); it != end; ++it) {
				w = *it;
				if (w == 0) {
					throw std::logic_error("elias_delta::encode(const SDSBitVector &v, SDSBitVector &z); entry of v equals 0 that cannot be encoded!");
				}
				// (number of bits to represent w) 
				len = bit_magic::l1BP(w) + 1;
				// (number of bits to represent the length of w) -1
				len_1_len = bit_magic::l1BP(len);
				// Write unary representation for the length of the length of w 
				bit_magic::writeIntAndMove(z_data, 1ULL << len_1_len, offset, len_1_len + 1);
				if (len_1_len) {
					bit_magic::writeIntAndMove(z_data, len, offset, len_1_len);
					bit_magic::writeIntAndMove(z_data, w, offset, len - 1);
				}
			}
			return true;
		}

		inline void elias_delta::encode(uint64_t x, uint64_t* &z, uint8_t &offset) {
			if (x == 0) {
				throw std::logic_error("elias_delta::encode(uint64_t x, uint64_t* &z, uint8_t &offset); x equals 0 that cannot be encoded!");
			}
			uint8_t len, len_1_len;
			// (number of bits to represent w) 
			len = bit_magic::l1BP(x) + 1;
			// (number of bits to represent the length of w) - 1
			len_1_len = bit_magic::l1BP(len);
			// Write unary representation for the length of the length of w 
			bit_magic::writeIntAndMove(z, 1ULL << len_1_len, offset, len_1_len + 1);
			if (len_1_len) {
				bit_magic::writeIntAndMove(z, len, offset, len_1_len);
				bit_magic::writeIntAndMove(z, x, offset, len - 1);
			}
		}

		template<class int_vector>
		bool elias_delta::decode(const int_vector &z, int_vector &v) {
			size_t len_1_len, len, n = 0;
			const uint64_t *z_data = z.data();
			const uint64_t *z_end = z.data() + (z.bit_size() >> 6);
			uint8_t offset = 0;
			while ((z_data < z_end) or (z_data == z_end and offset < (z.bit_size() & 0x3F))) {
				len_1_len = bit_magic::readUnaryIntAndMove(z_data, offset);
				//if(z_data==z.data() ) std::cerr<<" len_1_len="<<len_1_len<<" offset="<<(int)offset<<std::endl;		
				if (len_1_len) {
					len = bit_magic::readIntAndMove(z_data, offset, len_1_len) + (1 << len_1_len);
					//if(z_data==z.data() ) std::cerr<<" len="<<len<<" offset="<<(int)offset<<std::endl;		
					bit_magic::move_right(z_data, offset, len - 1);
				}
				++n;
			}
			v.setIntWidth(z.getIntWidth());
			//	n = z.m_elements;
			v.resize(n);
			return decode<false, true>(z.data(), 0, n, v.begin());
		}

		template<bool sumup, bool increment, class Iterator>
		inline uint64_t elias_delta::decode(const uint64_t *data, const size_type start_idx, size_type n, Iterator it) {
			data += (start_idx >> 6);
			uint64_t value = 0;
			size_type i = 0;
			size_type len_1_len, len;
			uint8_t offset = start_idx & 0x3F;
			while (i++ < n) {// while not all values are decoded
				len_1_len = bit_magic::readUnaryIntAndMove(data, offset); // read length of length of x
				if (!len_1_len) {
					value += 1;
				}
				else {
					len = bit_magic::readIntAndMove(data, offset, len_1_len) + (1ULL << len_1_len);
					value += bit_magic::readIntAndMove(data, offset, len - 1) + (1ULL << (len - 1));
				}
				if (increment) *(it++) = value;
				if (!sumup) value = 0;
			}
			return value;
		}
} // end namespace coder

#endif
