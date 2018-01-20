#ifndef _UTIL_HPP_
#define _UTIL_HPP_

#include "bitmagic.hpp"
#include <iosfwd> // forward declaration of ostream
#include <stdint.h> // for uint64_t uint32_t declaration
#include <iostream>// for cerr
#include <cassert>
#include <fstream> // file stream for storeToFile and loadFromFile
#include <ctime>  // for rand initialization
#include <string>

// macros to transform a defined name to a string
#define SDSL_STR(x) #x
#define SDSL_XSTR(s) SDSL_STR(s)	

	//! A namespace for helper functions
namespace util {
		//! Sets all bits of the int_vector to pseudo-random bits.
		template<class int_vector>
		void setRandomBits(int_vector &v);
		//! Sets all bits of the int_vector to 0-bits.
		template<class int_vector>
		void setZeroBits(int_vector &v);
		//! Sets all bits of the int_vector to 1-bits.
		template<class int_vector>
		void setOneBits(int_vector &v);
		//! Counts and returns the 1-bits an int_vector. 
		/*! \param v The int_vector to count the 1-bits.
		\return The number of 1-bits in v.
		*/
		template<class int_vector>
		typename int_vector::size_type getOneBits(const int_vector &v);

		//! Counts 10 bit pair occurencies. 
		/*! \sa getOneBits, getOneZeroBits
		*/
		template<class int_vector>
		typename int_vector::size_type getOneZeroBits(const int_vector &v);

		//! Counts 01 bit pair occurencies. 
		/*! \sa getOneBits, getZeroOneBits
		*/
		template <class int_vector>
		typename int_vector::size_type getZeroOneBits(const int_vector &v);

		//! Load a data structure from a file.
		/*! The data structure has to provide a load function.
		* \param v Data structure to load.
		\param file_name Name of the serialized file.
		*/
		template<class T>
		bool loadFromFile(T &v, const char *file_name);

		//! Store a data structure to a file.
		/*! The data structure has to provide a serialize function.
		\param v Data structure to store.
		\param file_name Name of the file where to store the data structure.
		*/
		template<class T>
		bool storeToFile(const T &v, const char *file_name);

		std::string demangle(const char* name);

		template<class T>
		typename T::size_type get_size_in_bytes(const T &v);

		struct nullstream : std::ostream {
			struct nullbuf : std::streambuf {
				int overflow(int c) { return traits_type::not_eof(c); }
			} m_sbuf;
			nullstream() : std::ios(&m_sbuf), std::ostream(&m_sbuf) {}
		};
}

	//==================== Template functions ====================

template<class T>
typename T::size_type util::get_size_in_bytes(const T &v) {
		if ((&v) == NULL)
			return 0;
		util::nullstream ns;
		return v.serialize(ns);
}

template<class T>
bool util::storeToFile(const T &v, const char *file_name) {
		std::ofstream out;
		out.open(file_name, std::ios::binary | std::ios::trunc | std::ios::out);
		if (!out)
			return false;
		v.serialize(out);
		out.close();
		return true;
}

template<class T>
bool util::loadFromFile(T &v, const char *file_name) {
		std::ifstream in;
		in.open(file_name, std::ios::binary | std::ios::in);
		if (!in)
			return false;
		v.load(in);
		in.close();
		return true;
}

template<class int_vector>
void util::setRandomBits(int_vector &v) {
		srand48((int)time(NULL));

		uint64_t *data = v.m_data;
		if (v.empty())
			return;
		*data = (((uint64_t)lrand48() & 0xFFFFULL) << 48)
			| (((uint64_t)lrand48() & 0xFFFFULL) << 32)
			| (((uint64_t)lrand48() & 0xFFFFULL) << 16)
			| ((uint64_t)lrand48() & 0xFFFFULL);
		for (typename int_vector::size_type i = 1; i < (v.capacity() >> 6); ++i) {
			*(++data) = (((uint64_t)lrand48() & 0xFFFFULL) << 48)
				| (((uint64_t)lrand48() & 0xFFFFULL) << 32)
				| (((uint64_t)lrand48() & 0xFFFFULL) << 16)
				| ((uint64_t)lrand48() & 0xFFFFULL);
		}
}

template<class int_vector>
void util::setZeroBits(int_vector &v) {
		uint64_t *data = v.m_data;
		if (v.empty())
			return;
		*data = 0ULL;
		for (typename int_vector::size_type i = 1; i < (v.capacity() >> 6); ++i) {
			*(++data) = 0ULL;
		}
}

template<class int_vector>
void util::setOneBits(int_vector &v) {
		uint64_t *data = v.m_data;
		if (v.empty())
			return;
		*data = 0xFFFFFFFFFFFFFFFFULL;
		for (typename int_vector::size_type i = 1; i < (v.capacity() >> 6); ++i) {
			*(++data) = 0xFFFFFFFFFFFFFFFFULL;
		}
}

template<class int_vector>
typename int_vector::size_type util::getOneBits(const int_vector &v) {
		const uint64_t *data = v.data();
		if (v.empty())
			return 0;
		typename int_vector::size_type result = (int_vector::size_type)bit_magic::b1Cnt(*data);
		for (typename int_vector::size_type i = 1; i < (v.capacity() >> 6); ++i) {
			result += (int_vector::size_type)bit_magic::b1Cnt(*(++data));
		}
		if (v.bit_size() & 0x3F) {
			result -= (int_vector::size_type)bit_magic::b1Cnt((*data) & (~bit_magic::Li1Mask[v.bit_size() & 0x3F]));
		}
		return result;
}

template<class int_vector>
typename int_vector::size_type util::getOneZeroBits(const int_vector &v) {
		const uint64_t *data = v.data();
		if (v.empty())
			return 0;
		uint64_t carry = 0, oldcarry = 0;
		typename int_vector::size_type result = bit_magic::b10Cnt(*data, carry);
		for (typename int_vector::size_type i = 1; i < (v.capacity() >> 6); ++i) {
			oldcarry = carry;
			result += (int_vector::size_type)bit_magic::b10Cnt(*(++data), carry);
		}
		if (v.bit_size() & 0x3F) {// if bit_size is not a multiple of 64, substract the counts of the additional bits
			result -= (int_vector::size_type)bit_magic::b1Cnt(bit_magic::b10Map(*data, oldcarry) & bit_magic::Li0Mask[v.bit_size() & 0x3F]);
		}
		return result;
}

template<class int_vector>
typename int_vector::size_type util::getZeroOneBits(const int_vector &v) {
		const uint64_t *data = v.data();
		if (v.empty())
			return 0;
		uint64_t carry = 1, oldcarry = 1;
		typename int_vector::size_type result = bit_magic::b01Cnt(*data, carry);
		for (typename int_vector::size_type i = 1; i < (v.capacity() >> 6); ++i) {
			oldcarry = carry;
			result += (int_vector::size_type)bit_magic::b01Cnt(*(++data), carry);
		}
		if (v.bit_size() & 0x3F) {// if bit_size is not a multiple of 64, substract the counts of the additional bits
			result -= (int_vector::size_type)bit_magic::b1Cnt(bit_magic::b01Map(*data, oldcarry) & bit_magic::Li0Mask[v.bit_size() & 0x3F]);
		}
		return result;
}

#endif // end file 

