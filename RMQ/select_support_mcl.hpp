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
/*! \file select_support_mcl.hpp
\brief select_support_mcl.hpp contains classes that support a sdsl::bit_vector with constant time select information.
\author Simon Gog
*/
#ifndef _SELECT_SUPPORT_JMC_HPP
#define _SELECT_SUPPORT_JMC_HPP

#include "vectors.hpp"
#include "select_support.hpp"

	template<uint8_t bit_pattern, uint8_t pattern_len>
	struct select_support_mcl_trait {
		typedef select_support::size_type	size_type;

		/* Count the number of arguments for the specific select support */
		static size_type arg_cnt(const bit_vector &v) {
			return 0;
		}

		static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry) {
			return 0;
		}

		static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry) {
			return 0;
		}

		static size_type args_in_the_word(uint64_t w, uint64_t &carry) {
			return 0;
		}

		static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry) {
			return 0;
		}

		static bool found_arg(size_type i, const bit_vector &v) {
			return 0;
		}

		static uint64_t init_carry(const uint64_t *data, size_type word_pos) {
			return 0;
		}
		static uint64_t get_carry(uint64_t w) {
			return 0;
		}
	};

	template<>
	struct select_support_mcl_trait<0, 1> {
		typedef select_support::size_type	size_type;

		static size_type arg_cnt(const bit_vector &v) {
			return v.bit_size() - util::getOneBits(v);
		}

		static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry) {
			return (size_type)bit_magic::b1Cnt((~w)& bit_magic::Li0Mask[offset]);
		}

		static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry) {
			return bit_magic::i1BP(~w & bit_magic::Li0Mask[offset], i);
		}
		static size_type args_in_the_word(uint64_t w, uint64_t &carry) {
			return (size_type)bit_magic::b1Cnt(~w);
		}
		static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry) {
			return bit_magic::i1BP(~w, i);
		}
		static bool found_arg(size_type i, const bit_vector &v) {
			return !v[i];
		}
		static uint64_t init_carry(const uint64_t *data, size_type word_pos) {
			return 0;
		}
		static uint64_t get_carry(uint64_t w) {
			return 0;
		}
	};

	template<>
	struct select_support_mcl_trait<1, 1> {
		typedef select_support::size_type	size_type;

		/* Count the number of arguments for the specific select support */
		static size_type arg_cnt(const bit_vector &v) {
			return util::getOneBits(v);
		}

		static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry) {
			return (size_type)bit_magic::b1Cnt(w & bit_magic::Li0Mask[offset]);
		}

		static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry) {
			return bit_magic::i1BP(w & bit_magic::Li0Mask[offset], i);
		}

		static size_type args_in_the_word(uint64_t w, uint64_t &carry) {
			return (size_type)bit_magic::b1Cnt(w);
		}

		static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry) {
			return bit_magic::i1BP(w, i);
		}

		static bool found_arg(size_type i, const bit_vector &v) {
			return v[i];
		}

		static uint64_t init_carry(const uint64_t *data, size_type word_pos) {
			return 0;
		}
		static uint64_t get_carry(uint64_t w) {
			return 0;
		}
	};

	template<>
	struct select_support_mcl_trait<10, 2> {
		typedef select_support::size_type	size_type;

		/* Count the number of arguments for the specific select support */
		static size_type arg_cnt(const bit_vector &v) {
			return util::getOneZeroBits(v);
		}

		static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry) {
			return (size_type)bit_magic::b1Cnt(bit_magic::b10Map(w, carry) & bit_magic::Li0Mask[offset]);
		}

		static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry) {
			return bit_magic::i1BP(bit_magic::b10Map(w, carry) & bit_magic::Li0Mask[offset], i);
		}

		static size_type args_in_the_word(uint64_t w, uint64_t &carry) {
			return bit_magic::b10Cnt(w, carry);
		}

		static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry) {
			return bit_magic::i1BP(bit_magic::b10Map(w, carry), i);
		}

		static bool found_arg(size_type i, const bit_vector &v) {
			if (i > 0 && v[i - 1] && !v[i])
				return true;
			return false;
		}

		static uint64_t init_carry(const uint64_t *data, size_type word_pos) {
			return word_pos ? (*(data - 1) >> 63) : 0;
		}
		static uint64_t get_carry(uint64_t w) {
			return w >> 63;
		}
	};

	template<>
	struct select_support_mcl_trait<01, 2> {
		typedef select_support::size_type	size_type;

		/* Count the number of arguments for the specific select support */
		static size_type arg_cnt(const bit_vector &v) {
			return util::getZeroOneBits(v);
		}

		static size_type args_in_the_first_word(uint64_t w, uint8_t offset, uint64_t carry) {
			return (size_type)bit_magic::b1Cnt(bit_magic::b01Map(w, carry) & bit_magic::Li0Mask[offset]);
		}

		static size_type ith_arg_pos_in_the_first_word(uint64_t w, size_type i, uint8_t offset, uint64_t carry) {
			return bit_magic::i1BP(bit_magic::b01Map(w, carry) & bit_magic::Li0Mask[offset], i);
		}

		static size_type args_in_the_word(uint64_t w, uint64_t &carry) {
			return bit_magic::b01Cnt(w, carry);
		}

		static size_type ith_arg_pos_in_the_word(uint64_t w, size_type i, uint64_t carry) {
			return bit_magic::i1BP(bit_magic::b01Map(w, carry), i);
		}

		static bool found_arg(size_type i, const bit_vector &v) {
			if (i > 0 && !v[i - 1] && v[i])
				return true;
			return false;
		}

		static uint64_t init_carry(const uint64_t *data, size_type word_pos) {
			return word_pos ? (*(data - 1) >> 63) : 1;
		}
		static uint64_t get_carry(uint64_t w) {
			return w >> 63;
		}
	};


	//! A class supporting constant time select queries (proposed by Munro/Clark, 1996) enhanced by broadword computing tricks.
	/*!
	* @ingroup select_support_group
	*/
	template<uint8_t b = 1, uint8_t pattern_len = 1>
	class select_support_mcl : public select_support {
	private:
		uint32_t m_logn, m_logn2, m_logn4; // log(Size of the Supported BitVector), logn2=(\log n)^2, logn4=(\log n)^4
		int_vector<0> m_superblock; // there exists \f$\frac{n}{4096}<\frac{n}{\log n}\f$ entries
									// entry i equals the answer to select_1(B,i*4096)
		int_vector<0> *m_longsuperblock;
		int_vector<0> *m_miniblock;
		size_type m_arg_cnt;
		void copy(const select_support_mcl<b, pattern_len> &ss);
		void construct();
		void initData();
	public:
		select_support_mcl(const int_vector<1> *v = NULL);
		select_support_mcl(const select_support_mcl<b, pattern_len> &ss);
		~select_support_mcl();
		void init(const int_vector<1> *v = NULL);
		void init_slow(const int_vector<1> *v = NULL);
		//! Select function
		/*! \sa select_support.select
		*/
		inline const size_type select(size_type i) const;
		//! Alias for select(i).
		inline const size_type operator()(size_type i)const;
		size_type serialize(std::ostream &out)const;
		void load(std::istream &in, const int_vector<1> *v = NULL);
		void set_vector(const int_vector<1> *v = NULL);
		select_support_mcl<b, pattern_len>& operator=(const select_support_mcl &ss);

		//! Swap operator
		/*! This swap Operator swaps two select_support_mcls in constant time.
		*  All members (excluded the pointer to the supported SDSBitVector) are swapped.
		*/
		void swap(select_support_mcl<b, pattern_len> &ss);
		//! Equality Operator
		/*! Two select_support_mcls are equal if all member variables are equal.
		* Required for the Equality Comparable Concept of the STL.
		* \sa operator!=
		*/
		bool operator==(const select_support_mcl<b, pattern_len> &ss)const;
		//! Unequality Operator
		/*! Two select_support_mcls are not equal if any member variable are not equal.
		* Required for the Equality Comparable Concept of the STL.
		* \sa operator==
		*/
		bool operator!=(const select_support_mcl<b, pattern_len> &ss)const;
	};


	template<uint8_t b, uint8_t pattern_len>
	void select_support_mcl<b, pattern_len>::construct() {
		m_longsuperblock = NULL;
		m_miniblock = NULL;
		m_arg_cnt = 0;
	}

	template<uint8_t b, uint8_t pattern_len>
	select_support_mcl<b, pattern_len>::select_support_mcl(const int_vector<1> *v) :select_support(v) {
		construct();
		init(v);
	}

	template<uint8_t b, uint8_t pattern_len>
	select_support_mcl<b, pattern_len>::select_support_mcl(const select_support_mcl &ss) :select_support(ss.m_v) {
		construct();
		copy(ss);
	}

	template<uint8_t b, uint8_t pattern_len>
	select_support_mcl<b, pattern_len>& select_support_mcl<b, pattern_len>::operator=(const select_support_mcl &ss) {
		if (this != &ss) {
			copy(ss);
		}
		return *this;
	}

	template<uint8_t b, uint8_t pattern_len>
	void select_support_mcl<b, pattern_len>::swap(select_support_mcl &ss) {
		std::swap(m_logn, ss.m_logn);
		std::swap(m_logn2, ss.m_logn2);
		std::swap(m_logn4, ss.m_logn4);
		m_superblock.swap(ss.m_superblock);
		std::swap(m_longsuperblock, ss.m_longsuperblock);
		std::swap(m_miniblock, ss.m_miniblock);
		std::swap(m_arg_cnt, ss.m_arg_cnt);
		//	std::swap(m_v, ss.m_v);
	}

	template<uint8_t b, uint8_t pattern_len>
	void select_support_mcl<b, pattern_len>::copy(const select_support_mcl<b, pattern_len> &ss) {
		m_logn = ss.m_logn;		// copy log n
		m_logn2 = ss.m_logn2;		// copy (logn)^2
		m_logn4 = ss.m_logn4;		// copy (logn)^4
		m_superblock = ss.m_superblock;  // copy long superblock
		m_arg_cnt = ss.m_arg_cnt;	    // copy count of 1-bits
		m_v = ss.m_v;		    // copy pointer to the supported bit vector
		size_type sb = (m_arg_cnt + 4095) >> 12;
		if (m_longsuperblock != NULL)
			delete[] m_longsuperblock;
		m_longsuperblock = NULL;
		if (ss.m_longsuperblock != NULL) {
			m_longsuperblock = new int_vector<0>[sb]; //copy longsuperblocks
			for (size_type i = 0; i<sb; ++i) {
				m_longsuperblock[i] = ss.m_longsuperblock[i];
			}
		}
		if (m_miniblock != NULL)
			delete[] m_miniblock;
		m_miniblock = NULL;
		if (ss.m_miniblock != NULL) {
			m_miniblock = new int_vector<0>[sb]; // copy miniblocks
			for (size_type i = 0; i<sb; ++i) {
				m_miniblock[i] = ss.m_miniblock[i];
			}
		}
	}

	template<uint8_t b, uint8_t pattern_len>
	select_support_mcl<b, pattern_len>::~select_support_mcl() {
		if (m_longsuperblock != NULL)
			delete[] m_longsuperblock;
		if (m_miniblock != NULL)
			delete[] m_miniblock;
	}

	template<uint8_t b, uint8_t pattern_len>
	void select_support_mcl<b, pattern_len>::init_slow(const int_vector<1> *v) {
		set_vector(v);
		initData();
		if (m_v == NULL)
			return;
		// Count the number of arguments in the bit vector
		m_arg_cnt = select_support_mcl_trait<b, pattern_len>::arg_cnt(*v);

		const size_type SUPER_BLOCK_SIZE = 4096;

		if (m_arg_cnt == 0) // if there are no arguments in the vector we are done...
			return;

		size_type sb = (m_arg_cnt + SUPER_BLOCK_SIZE - 1) / SUPER_BLOCK_SIZE; // number of superblocks
		if (m_miniblock != NULL) delete[] m_miniblock;
		m_miniblock = new int_vector<0>[sb];

		m_superblock = int_vector<0>(sb, 0, m_logn);// TODO: hier koennte man logn noch optimieren...s


		bit_vector::size_type arg_position[SUPER_BLOCK_SIZE], arg_cnt = 0;
		for (size_type i = 0, sb_cnt = 0; i < v->size(); ++i) {
			if (select_support_mcl_trait<b, pattern_len>::found_arg(i, *v)) {
				arg_position[arg_cnt%SUPER_BLOCK_SIZE] = i;
				assert(arg_position[arg_cnt%SUPER_BLOCK_SIZE] == i);
				++arg_cnt;
				if (arg_cnt % SUPER_BLOCK_SIZE == 0 || arg_cnt == m_arg_cnt) { // 
					assert(sb_cnt < sb);
					m_superblock[sb_cnt] = arg_position[0];


					size_type pos_diff = arg_position[(arg_cnt - 1) % SUPER_BLOCK_SIZE] - arg_position[0];
					if (pos_diff > m_logn4) {// longblock
						if (m_longsuperblock == NULL) m_longsuperblock = new int_vector<0>[sb]; // create longsuperblock
						m_longsuperblock[sb_cnt] = int_vector<0>(SUPER_BLOCK_SIZE, 0, bit_magic::l1BP(arg_position[(arg_cnt - 1) % SUPER_BLOCK_SIZE]) + 1);

						for (size_type j = 0; j <= (arg_cnt - 1) % SUPER_BLOCK_SIZE; ++j) m_longsuperblock[sb_cnt][j] = arg_position[j]; // copy argument positions to longsuperblock
					}
					else {// short block
						m_miniblock[sb_cnt] = int_vector<0>(64, 0, bit_magic::l1BP(pos_diff) + 1);
						for (size_type j = 0; j <= (arg_cnt - 1) % SUPER_BLOCK_SIZE; j += 64) {
							m_miniblock[sb_cnt][j / 64] = arg_position[j] - arg_position[0];
						}
					}
					++sb_cnt;
				}
			}
		}

	}

	template<uint8_t b, uint8_t pattern_len>
	void select_support_mcl<b, pattern_len>::init(const int_vector<1> *v) {
		init_slow(v);
		return;
	}

	template<uint8_t b, uint8_t pattern_len>
	inline const typename select_support_mcl<b, pattern_len>::size_type select_support_mcl<b, pattern_len>::select(size_type i)const {
		assert(i > 0 && i <= m_arg_cnt);

		i = i - 1;
		size_type sb_idx = i >> 12;   // i/4096
		size_type offset = i & 0xFFF; // i%4096
		if (m_longsuperblock != NULL && !m_longsuperblock[sb_idx].empty()) {
			return (select_support_mcl<b, pattern_len>::size_type)m_longsuperblock[sb_idx][offset];
		}
		else {
			if ((offset & 0x3F) == 0) {
				assert(sb_idx < m_superblock.size());
				if ((offset >> 6) >= m_miniblock[sb_idx].size()) {
					std::cerr << " i=" << i << std::endl;
					std::cerr << " " << (offset >> 6) << " >= " << m_miniblock[sb_idx].size() << std::endl;
				}
				assert((offset >> 6) < m_miniblock[sb_idx].size());
				return (select_support_mcl<b, pattern_len>::size_type)(m_superblock[sb_idx] + m_miniblock[sb_idx][offset >> 6/*/64*/]);
			}
			else {
				i = i - (sb_idx << 12) - ((offset >> 6) << 6);
				// now i > 0 and i <= 64
				assert(i > 0);

				size_type pos = (size_type)(m_superblock[sb_idx] + m_miniblock[sb_idx][offset >> 6] + 1);

				// now pos is the position from where we search for the ith argument 
				size_type word_pos = pos >> 6;
				size_type word_off = pos & 0x3F;

				const uint64_t *data = m_v->data() + word_pos;
				uint64_t carry = select_support_mcl_trait<b, pattern_len>::init_carry(data, word_pos);
				size_type args = select_support_mcl_trait<b, pattern_len>::args_in_the_first_word(*data, (uint8_t)word_off, carry);

				if (args >= i) {
					return (word_pos << 6) + select_support_mcl_trait<b, pattern_len>::ith_arg_pos_in_the_first_word(*data, i, (uint8_t)word_off, carry);
				}
				word_pos += 1;
				size_type sum_args = args;
				carry = select_support_mcl_trait<b, pattern_len>::get_carry(*data);
				uint64_t old_carry = carry;
				args = select_support_mcl_trait<b, pattern_len>::args_in_the_word(*(++data), carry);
				while (sum_args + args < i) {
					sum_args += args;
					assert(data + 1 < m_v->data() + (m_v->capacity() >> 6));
					old_carry = carry;
					args = select_support_mcl_trait<b, pattern_len>::args_in_the_word(*(++data), carry);
					word_pos += 1;
				}
				return (word_pos << 6) + select_support_mcl_trait<b, pattern_len>::ith_arg_pos_in_the_word(*data, i - sum_args, old_carry);
			}
		}
	}

	template<uint8_t b, uint8_t pattern_len>
	inline const typename select_support_mcl<b, pattern_len>::size_type select_support_mcl<b, pattern_len>::operator()(size_type i)const {
		return select(i);
	}

	template<uint8_t b, uint8_t pattern_len>
	void select_support_mcl<b, pattern_len>::initData() {
		m_arg_cnt = 0;
		if (m_v == NULL) {
			m_logn = m_logn2 = m_logn4 = 0;
		}
		else {
			m_logn = bit_magic::l1BP(m_v->capacity()) + 1;
			m_logn2 = m_logn*m_logn;
			m_logn4 = m_logn2*m_logn2;
		}
		if (m_longsuperblock != NULL)
			delete[] m_longsuperblock;
		m_longsuperblock = NULL;
		if (m_miniblock != NULL)
			delete[] m_miniblock;
		m_miniblock = NULL;
	}

	template<uint8_t b, uint8_t pattern_len>
	void select_support_mcl<b, pattern_len>::set_vector(const int_vector<1> *v) {
		m_v = v;
	}

	template<uint8_t b, uint8_t pattern_len>
	typename select_support_mcl<b, pattern_len>::size_type select_support_mcl<b, pattern_len>::serialize(std::ostream &out)const {
		size_type written_bytes = 0;
		// write the number of 1-bits in the supported SDSBitVector
		out.write((char *)&m_arg_cnt, sizeof(size_type) / sizeof(char));
		written_bytes = sizeof(size_type) / sizeof(char);
		// number of superblocks in the data structure
		size_type sb = (m_arg_cnt + 4095) >> 12;

		if (m_arg_cnt) { // if there exists 1-bits to be supported
			written_bytes += m_superblock.serialize(out); // serialize superblocks
			int_vector<1> mini_or_long;// Helper vector: mini or long block?
			if (m_longsuperblock != NULL) {
				mini_or_long.resize(sb); // resize indicatior BitVector to the number of superblocks
				for (size_type i = 0; i< sb; ++i)
					mini_or_long[i] = !m_miniblock[i].empty();
			}
			written_bytes += mini_or_long.serialize(out);
			for (size_type i = 0; i < sb; ++i)
				if (!mini_or_long.empty() && !mini_or_long[i]) {
					written_bytes += m_longsuperblock[i].serialize(out);
					//		   			SDSCoder::compress(m_longsuperblock[i], z);
					//					z.serialize(out);
					//SDSBitVector *z = SDSCoder::compress(m_longsuperblock[i]);
					//if(z!=NULL) delete z;
				}
				else {
					written_bytes += m_miniblock[i].serialize(out);
					//		   			SDSCoder::compress(m_miniblock[i], z);
					//		   			SDSCoder::decompress(z, zz);
					//					assert(m_miniblock[i]==zz);
					//					z.serialize(out);
					//SDSBitVector *z = SDSCoder::compress(m_miniblock[i]);
					//if(z!=NULL) delete z;
				}
		}
		return written_bytes;
	}

	template<uint8_t b, uint8_t pattern_len>
	void select_support_mcl<b, pattern_len>::load(std::istream &in, const int_vector<1> *v) {
		set_vector(v);
		initData();
		// read the number of 1-bits in the supported SDSBitVector
		in.read((char *)&m_arg_cnt, sizeof(size_type) / sizeof(char));
		size_type sb = (m_arg_cnt + 4095) >> 12;

		if (m_arg_cnt) { // if there exists 1-bits to be supported
						 //			SDSBitVector z;
						 //			z.load(in);// load compressed superblocks
						 //			SDSCoder::decompress(z, m_superblock);
			m_superblock.load(in); // load superblocks

			if (m_miniblock != NULL) {
				delete[] m_miniblock;
				m_miniblock = NULL;
			}
			if (m_longsuperblock != NULL) {
				delete[] m_longsuperblock;
				m_longsuperblock = NULL;
			}

			int_vector<1> mini_or_long;// Helper vector: mini or long block?
			mini_or_long.load(in); // Load the helper vector
			m_miniblock = new int_vector<0>[sb]; // Create miniblock int_vector<0>
			if (!mini_or_long.empty())
				m_longsuperblock = new int_vector<0>[sb]; // Create longsuperblock int_vector<0>

			for (size_type i = 0; i < sb; ++i)
				if (!mini_or_long.empty() && !mini_or_long[i]) {
					m_longsuperblock[i].load(in);
					//					z.load(in);
					//					SDSCoder::decompress(z,m_longsuperblock[i]);
				}
				else {
					m_miniblock[i].load(in);
					//					z.load(in);
					//					SDSCoder::decompress(z,m_miniblock[i]);
				}
		}
	}

	template<uint8_t b, uint8_t pattern_len>
	bool select_support_mcl<b, pattern_len>::operator==(const select_support_mcl &ss)const {
		if (this == &ss)
			return true;
		if (m_logn == ss.m_logn and m_logn2 == ss.m_logn2 and m_logn4 == ss.m_logn4
			and m_arg_cnt == ss.m_arg_cnt and  m_superblock == ss.m_superblock) {
			if (m_longsuperblock == NULL and ss.m_longsuperblock != NULL)
				return false;
			if (ss.m_longsuperblock == NULL and m_longsuperblock != NULL)
				return false;
			if (m_miniblock == NULL and ss.m_miniblock != NULL)
				return false;
			if (ss.m_miniblock == NULL and m_miniblock != NULL)
				return false;
			size_type sb = (m_arg_cnt + 4095) >> 12;
			if (m_miniblock != NULL)
				for (size_type i = 0; i < sb; ++i)
					if (m_miniblock[i] != ss.m_miniblock[i])
						return false;
			if (m_longsuperblock != NULL)
				for (size_type i = 0; i < sb; ++i)
					if (m_longsuperblock[i] != ss.m_longsuperblock[i])
						return false;
			if (m_v == NULL and ss.m_v != NULL)
				return false;
			if (m_v != NULL and ss.m_v == NULL)
				return false;
			if (m_v == NULL and ss.m_v == NULL)
				return false;
			return *ss.m_v == *m_v;// supported vectors have to be equal. See docu!!!
		}
		else
			return false;
	}

	template<uint8_t b, uint8_t pattern_len>
	bool select_support_mcl<b, pattern_len>::operator!=(const select_support_mcl &ss)const {
		return !(*this == ss);
	}

#endif
