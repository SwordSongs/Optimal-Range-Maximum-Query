/* sdsl - succinct data structures library
Copyright (C) 2008 Simon Gog

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
/*! \file rank_support_v.hpp
\brief rank_support_v.hpp contains rank_support_v that support a sdsl::bit_vector with constant time rank information.
\author Simon Gog
*/
#ifndef _RANK_SUPPORT_V
#define RANK_SUPPORT_V

#include "rank_support.hpp"
#include "vectors.hpp"

template<uint8_t bit_pattern, uint8_t pattern_len>
struct rank_support_v_trait {
		typedef rank_support::size_type	size_type;

		static size_type args_in_the_word(uint64_t w, uint64_t &carry) {
			return 0;
		}

		static uint32_t word_rank(const uint64_t *data, size_type idx) {
			return 0;
		}

		static uint64_t init_carry() {
			return 0;
		}
};

template<>
struct rank_support_v_trait<0, 1> {
		typedef rank_support::size_type	size_type;

		static size_type args_in_the_word(uint64_t w, uint64_t &carry) {
			return (size_type)bit_magic::b1Cnt(~w);
		}

		static uint32_t word_rank(const uint64_t* data, size_type idx) {
			return	(uint32_t)bit_magic::b1Cnt((~*(data + (idx >> 6))) & bit_magic::Li1Mask[idx & 0x3F]);
		}

		static uint64_t init_carry() {
			return 0;
		}
};

template<>
struct rank_support_v_trait<1, 1> {
		typedef rank_support::size_type	size_type;

		static size_type args_in_the_word(uint64_t w, uint64_t &carry) {
			return (size_type)bit_magic::b1Cnt(w);
		}

		static uint32_t word_rank(const uint64_t* data, size_type idx) {
			return	(uint32_t)bit_magic::b1Cnt(*(data + (idx >> 6)) & bit_magic::Li1Mask[idx & 0x3F]);
		}

		static uint64_t init_carry() {
			return 0;
		}
};

template<>
struct rank_support_v_trait<10, 2> {
		typedef rank_support::size_type	size_type;

		static size_type args_in_the_word(uint64_t w, uint64_t &carry) {
			return bit_magic::b10Cnt(w, carry);
		}

		static uint32_t word_rank(const uint64_t* data, size_type idx) {
			data = data + (idx >> 6);
			uint64_t carry = (idx>63) ? *(data - 1) >> 63 : 0;
			return	(uint32_t)bit_magic::b1Cnt(bit_magic::b10Map(*data, carry) & bit_magic::Li1Mask[idx & 0x3F]);
		}

		static uint64_t init_carry() {
			return 0;
		}
};

template<>
struct rank_support_v_trait<01, 2> {
		typedef rank_support::size_type	size_type;

		static size_type args_in_the_word(uint64_t w, uint64_t &carry) {
			return bit_magic::b01Cnt(w, carry);
		}

		static uint32_t word_rank(const uint64_t* data, size_type idx) {
			data = data + (idx >> 6);
			uint64_t carry = (idx>63) ? *(data - 1) >> 63 : 0;
			return	(uint32_t)bit_magic::b1Cnt(bit_magic::b01Map(*data, carry) & bit_magic::Li1Mask[idx & 0x3F]);
		}

		static uint64_t init_carry() {
			return 1;
		}
};

	//! A class supporting rank queries in constant time. The implementation is a version of the data structure porposed by Vigna (WEA 2008).
	/*! \par Space complexity
	*  \f$ 0.25n\f$ for a bit vector of length n bits.
	* @ingroup rank_support_group
	*/
template<uint8_t b = 1, uint8_t pattern_len = 1>
class rank_support_v : public rank_support {
	private:
		int_vector<64> m_basic_block; // basic block for interleaved storage of superblockrank and blockrank
	public:
		rank_support_v(const int_vector<1> *v = NULL);
		rank_support_v(const rank_support_v<b, pattern_len> &rs);
		~rank_support_v();
		void init(const int_vector<1> *v = NULL);
		const size_type rank(size_type idx) const;
		const size_type operator()(size_type idx)const;
		size_type serialize(std::ostream &out)const;
		void load(std::istream &in, const int_vector<1> *v = NULL);
		void set_vector(const int_vector<1> *v = NULL);

		//! Assign Operator
		/*! Required for the Assignable Concept of the STL.
		*/
		rank_support_v& operator=(const rank_support_v<b, pattern_len> &rs);
		//! swap Operator
		/*! Swap two rank_support_v in constant time.
		*	All members (excluded the pointer to the supported SDSBitVector) are swapped.
		*
		*  Required for the Container Concept of the STL.
		*/
		void swap(rank_support_v<b, pattern_len> &rs);
		//! Equality Operator
		/*! Two rank_support_vs are equal if all member variables are equal.
		*
		* Required for the Equality Comparable Concept of the STL.
		* \sa operator!=
		*/
		bool operator==(const rank_support_v<b, pattern_len> &rs)const;
		//! Unequality Operator
		/*! Two rank_support_vs are not equal if any member variable are not equal.
		*
		* Required for the Equality Comparable Concept of the STL.
		* \sa operator==
		*/
		bool operator!=(const rank_support_v<b, pattern_len> &rs)const;

};

template<uint8_t b, uint8_t pattern_len>
inline rank_support_v<b, pattern_len>::rank_support_v(const int_vector<1> *v) {
		init(v);
}

template<uint8_t b, uint8_t pattern_len>
inline rank_support_v<b, pattern_len>::rank_support_v(const rank_support_v &rs) {
		m_v = rs.m_v;
		m_basic_block = rs.m_basic_block;
}

template<uint8_t b, uint8_t pattern_len>
inline void rank_support_v<b, pattern_len>::init(const int_vector<1> *v) {
		set_vector(v);
		if (v == NULL || v->empty())
			return;
		size_type basic_block_size = ((v->capacity() >> 9) + 1) << 1;
		m_basic_block.resize(basic_block_size); // resize structure for basic_blocks
		if (m_basic_block.empty())
			return;
		const uint64_t *data = m_v->data();
		size_type i, j = 0;
		m_basic_block[0] = m_basic_block[1] = 0;
		//	uint64_t carry = 0;

		uint64_t carry = rank_support_v_trait<b, pattern_len>::init_carry();
		uint64_t sum = rank_support_v_trait<b, pattern_len>::args_in_the_word(*data, carry)/*bit_magic::b1Cnt(*data)*/, second_level_cnt = 0;
		for (i = 1; i < (m_v->capacity() >> 6); ++i) {
			if (!(i & 0x7)) {// if i%8==0
				j += 2;
				m_basic_block[j - 1] = second_level_cnt;
				m_basic_block[j] = m_basic_block[j - 2] + sum;
				second_level_cnt = sum = 0;
			}
			else
			{
				second_level_cnt |= sum << (63 - 9 * (i & 0x7));//  54, 45, 36, 27, 18, 9, 0
			}
			sum += rank_support_v_trait<b, pattern_len>::args_in_the_word(*(++data), carry);
			//		sum += bit_magic::b1Cnt(*(++data));
		}
		if (i & 0x7) { // if i%8 != 0
			second_level_cnt |= sum << (63 - 9 * (i & 0x7));
			m_basic_block[j + 1] = second_level_cnt;
		}
		else { // if i%8 == 0
			j += 2;
			m_basic_block[j - 1] = second_level_cnt;
			m_basic_block[j] = m_basic_block[j - 2] + sum;
			m_basic_block[j + 1] = 0;
		}
}

template<uint8_t b, uint8_t pattern_len>
inline const typename rank_support_v<b, pattern_len>::size_type rank_support_v<b, pattern_len>::rank(size_type idx)const {
		const uint64_t *p = m_basic_block.data() + ((idx >> 8) & 0xFFFFFFFFFFFFFFFEULL);// (idx/512)*2
		if (idx & 0x3F)
			return  (const rank_support::size_type)(*p + ((*(p + 1) >> (63 - 9 * ((idx & 0x1FF) >> 6))) & 0x1FF) +
			rank_support_v_trait<b, pattern_len>::word_rank(m_v->data(), idx));
		//		   		bit_magic::b1Cnt( (*(m_v->data()+(idx>>6)))&bit_magic::Li1Mask[idx&0x3F] );
		else
			return  (const rank_support::size_type)(*p + ((*(p + 1) >> (63 - 9 * ((idx & 0x1FF) >> 6))) & 0x1FF));
}


template<uint8_t b, uint8_t pattern_len>
inline const typename rank_support_v<b, pattern_len>::size_type rank_support_v<b, pattern_len>::operator()(size_type idx)const {
		return rank(idx);
}

template<uint8_t b, uint8_t pattern_len>
inline rank_support_v<b, pattern_len>::~rank_support_v() {}

template<uint8_t b, uint8_t pattern_len>
inline void rank_support_v<b, pattern_len>::set_vector(const int_vector<1> *v) {
		m_v = v;
}

template<uint8_t b, uint8_t pattern_len>
inline typename rank_support_v<b, pattern_len>::size_type rank_support_v<b, pattern_len>::serialize(std::ostream &out)const {
		return m_basic_block.serialize(out);
}

template<uint8_t b, uint8_t pattern_len>
inline void rank_support_v<b, pattern_len>::load(std::istream &in, const int_vector<1> *v) {
		set_vector(v);
		assert(m_v != NULL); // supported bit vector should be known
		m_basic_block.load(in);
}

template<uint8_t b, uint8_t pattern_len>
inline rank_support_v<b, pattern_len>& rank_support_v<b, pattern_len>::operator=(const rank_support_v &rs) {
		if (this != &rs) {
			set_vector(rs.m_v);
			m_basic_block = rs.m_basic_block;
		}
		return *this;
}

template<uint8_t b, uint8_t pattern_len>
inline void rank_support_v<b, pattern_len>::swap(rank_support_v &rs) {
		if (this != &rs) { // if rs and _this_ are not the same object
						   // TODO: swap m_v??? no!!!
						   //		std::swap(m_v, rs.m_v);
			m_basic_block.swap(rs.m_basic_block);
		}
}

template<uint8_t b, uint8_t pattern_len>
inline bool rank_support_v<b, pattern_len>::operator==(const rank_support_v &rs)const {
		if (this == &rs)
			return true;
		return m_basic_block == rs.m_basic_block and *(rs.m_v) == *m_v;
}

template<uint8_t b, uint8_t pattern_len>
inline bool rank_support_v<b, pattern_len>::operator!=(const rank_support_v &rs)const {
		return !(*this == rs);
}

#endif // end file 
