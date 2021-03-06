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
/*! \file int_vector.hpp
\brief int_vector.hpp contains the sdsl::int_vector class.
\author Simon Gog
*/
#ifndef _VECTORS_HPP_
#define _VECTORS_HPP_

#include "bitmagic.hpp"
#include "util.hpp"
#include <iosfwd> // forward declaration of ostream
#include <stdint.h> // for uint64_t uint32_t declaration
#include <stdexcept> // for exceptions
#include <iostream>// for cerr
#include <cassert>
#include <iterator>
#include <cstdlib>
#include <ctime>  // for rand initialization
#include <sstream>
#include <cstring> // for memcpy
#include <ostream>
#include <istream>
#include <string>

template<uint8_t fixedIntWidth = 0, class size_type_class = size_t>
class int_vector; // forward declaration

					  //! bit_vector is a specialization of the int_vector. 
typedef int_vector<1> bit_vector;

template<class int_vector>
class int_vector_reference;

template<class int_vector>
struct int_vector_iterator_base;	 // forward declaration

									 //class SDSBitVectorIterator;  // forward declaration
template<class int_vector>
struct int_vector_iterator; // forward declaration

template<class int_vector>
struct int_vector_const_iterator;  // forward declaration

class rank_support_jmc;
template<uint8_t bit_pattern, uint8_t pattern_len>
class rank_support_v;
class rank_support;

class select_support;
template<class RankSupport>
class select_support_bs;
template<uint8_t bit_pattern, uint8_t pattern_len>
class select_support_mcl;

namespace coder {
		class fibonacci;
		class elias_delta;
		class ternary;
}


template<uint8_t fixedIntWidth, class size_type_class>
struct int_vector_trait {
		typedef uint64_t   				value_type;
		typedef size_type_class			size_type;
		typedef int_vector_reference<int_vector<fixedIntWidth> > 	reference;
		typedef const uint64_t const_reference;
		typedef uint8_t int_width_type;
		typedef int_vector_iterator<int_vector<fixedIntWidth> > 		iterator;
		typedef int_vector_const_iterator<int_vector<fixedIntWidth> >   const_iterator;
		// Sets int_width to new_int_width
		static void setIntWidth(int_width_type &int_width, uint8_t new_int_width) {
			if (fixedIntWidth == 1)
				return;
			if (new_int_width>0 && new_int_width <= 64)
				int_width = new_int_width;
			else
				int_width = 64;
		}

		static iterator begin(int_vector<fixedIntWidth> *v, uint64_t *begin) {
			return iterator(v, 0);
		}
		static iterator end(int_vector<fixedIntWidth> *v, uint64_t *begin, size_type size) {
			return iterator(v, v->size()*v->getIntWidth());
		}
		static const_iterator begin(const int_vector<fixedIntWidth> *v, const uint64_t *begin) {
			return const_iterator(v, 0);
		}
		static const_iterator end(const int_vector<fixedIntWidth> *v, const uint64_t *begin, size_type size) {
			return const_iterator(v, v->size()*v->getIntWidth());
		}
};

template<class size_type_class>
struct int_vector_trait<64, size_type_class> {
		typedef uint64_t   				value_type;
		typedef size_type_class			size_type;
		typedef uint64_t& reference;
		typedef const uint64_t const_reference;
		typedef const uint8_t int_width_type;
		typedef uint64_t* iterator;
		typedef const uint64_t* const_iterator;


		static void setIntWidth(int_width_type &int_width, uint8_t new_int_width) {}

		static iterator begin(int_vector<64> *v, uint64_t *begin) {
			return begin;
		}
		static iterator end(int_vector<64> *v, uint64_t *begin, size_type size) {
			return begin + size;
		}
		static const_iterator begin(const int_vector<64> *v, const uint64_t *begin) {
			return begin;
		}
		static const_iterator end(const int_vector<64> *v, const uint64_t *begin, size_type size) {
			return begin + size;
		}
};

template<class size_type_class>
struct int_vector_trait<32, size_type_class> {
		typedef uint32_t   				value_type;
		typedef size_type_class			size_type;
		typedef uint32_t& reference;
		typedef const uint32_t const_reference;
		typedef const uint8_t int_width_type;
		typedef uint32_t* iterator;
		typedef const uint32_t* const_iterator;
		static void setIntWidth(int_width_type &int_width, uint8_t new_int_width) {}

		static iterator begin(int_vector<32> *v, uint64_t *begin) {
			return (uint32_t*)begin;
		}
		static iterator end(int_vector<32> *v, uint64_t *begin, size_type size) {
			return ((uint32_t*)begin) + size;
		}
		static const_iterator begin(const int_vector<32> *v, const uint64_t *begin) {
			return (uint32_t*)begin;
		}
		static const_iterator end(const int_vector<32> *v, const uint64_t *begin, size_type size) {
			return ((uint32_t*)begin) + size;
		}
};

template<class size_type_class>
struct int_vector_trait<16, size_type_class> {
		typedef uint16_t   				value_type;
		typedef size_type_class			size_type;
		typedef uint16_t& reference;
		typedef const uint16_t const_reference;
		typedef const uint8_t int_width_type;
		typedef uint16_t* iterator;
		typedef const uint16_t* const_iterator;
		static void setIntWidth(int_width_type &int_width, uint8_t new_int_width) {}

		static iterator begin(int_vector<16> *v, uint64_t *begin) {
			return (uint16_t*)begin;
		}
		static iterator end(int_vector<16> *v, uint64_t *begin, size_type size) {
			return ((uint16_t*)begin) + size;
		}
		static const_iterator begin(const int_vector<16> *v, const uint64_t *begin) {
			return (uint16_t*)begin;
		}
		static const_iterator end(const int_vector<16> *v, const uint64_t *begin, size_type size) {
			return ((uint16_t*)begin) + size;
		}
};

template<class size_type_class>
struct int_vector_trait<8, size_type_class> {
		typedef uint8_t					value_type;
		typedef size_type_class			size_type;
		typedef uint8_t& reference;
		typedef const uint8_t const_reference;
		typedef const uint8_t int_width_type;
		typedef uint8_t* iterator;
		typedef const uint8_t* const_iterator;
		static void setIntWidth(int_width_type &int_width, uint8_t new_int_width) {}

		static iterator begin(int_vector<8> *v, uint64_t *begin) {
			return (uint8_t*)begin;
		}
		static iterator end(int_vector<8> *v, uint64_t *begin, size_type size) {
			return ((uint8_t*)begin) + size;
		}
		static const_iterator begin(const int_vector<8> *v, const uint64_t *begin) {
			return (uint8_t*)begin;
		}
		static const_iterator end(const int_vector<8> *v, const uint64_t *begin, size_type size) {
			return ((uint8_t*)begin) + size;
		}
};

	//! A generic vector class for integers of width \f$w\in [1..64]\f$.
	/*! \author Simon Gog

	This generic vector class could be used to generate a vector
	that contains integers of fixed width \f$w\in [1..64]\f$.
	E.g.
	\code
	sdsl::int_vector<34> v(10);
	\endcode
	creates a vector that can hold ten integers of width 34. You
	could use the object like a container of the Standard Template Library
	(STL, see http://www.sgi.com/tech/stl/ )
	as we have implemented all requirements for this concept. Therefore
	you could apply many algorithms of the stl to our datastructure.
	Let's have a look at bit_vector_and_stl.e.cpp for an example.

	If you set the fixed width parameter to zero
	\code
	sdsl::int_vector<0> v(10);
	\endcode
	you get a vector v of 10 elements which can each hold a 64-bit integer by default.
	See the int_vector_variable.e.cpp example for the behavior of the vector if
	you change the int width by calling the setIntWidth() method.

	*/
template<uint8_t fixedIntWidth, class size_type_class>
class int_vector {
	public:

		typedef typename int_vector_trait<fixedIntWidth, size_type_class>::value_type		value_type;  	// STL Container requirement
		typedef typename int_vector_trait<fixedIntWidth, size_type_class>::iterator 			iterator;    	// STL Container requirement
		typedef typename int_vector_trait<fixedIntWidth, size_type_class>::const_iterator	const_iterator;
		typedef typename int_vector_trait<fixedIntWidth, size_type_class>::reference 		reference;
		typedef typename int_vector_trait<fixedIntWidth, size_type_class>::const_reference	const_reference;
		typedef int_vector_reference<int_vector>*							pointer; // TODO: fix this for 63,32,16,8 bit vector
		typedef const value_type*											const_pointer;
		typedef ptrdiff_t 					difference_type;// STL Container requirement
		typedef size_type_class												size_type;		// STL Container requirement

		friend struct int_vector_trait<fixedIntWidth, size_type_class>;
		friend struct int_vector_iterator_base<int_vector>;
		friend struct int_vector_iterator<int_vector>;
		friend struct int_vector_const_iterator<int_vector>;
		friend class  coder::elias_delta;
		friend class  coder::fibonacci;
		friend class  coder::ternary;
		//! Operator to create an int_vector<1> (aka bit_vector) from an input stream.
		friend std::istream& operator>>(std::istream&, int_vector<1> &);
		friend void util::setRandomBits<int_vector>(int_vector &v);
		friend void util::setZeroBits<int_vector>(int_vector &v);
		friend void util::setOneBits<int_vector>(int_vector &v);
	private:
		size_type	m_size; //!< Number of bits needed to store int_vector.
		uint64_t *  m_data; //!< Pointer to the memory for the bits.
		typename int_vector_trait<fixedIntWidth, size_type_class>::int_width_type m_intWidth;//!< Width of the integers that are accessed via the [] operator .
	public:

		//! Constructor for int_vector.
		/*! \param elements The number of elements in the int_vector. Default value is 0.
		\param default_value The default value to initialize the elements.
		\param intWidth The width of integers which could be accessed via the [] operator.
		\sa resize, setIntWidth
		*/
		int_vector(size_type elements = 0, value_type default_value = 0, uint8_t intWidth = fixedIntWidth);

		//! Copy constructor for int_vector.
		/*! \param v The int_vector to copy
		Required for the STL Assignable Concept
		*/
		int_vector(const int_vector &v);

		//! Destructor for int_vector.
		~int_vector();

		//!	Equivalent to size() == 0.
		/*! Required for the STL Container Concept
		\sa size()
		*/
		bool empty() const;

		//! Swap method for int_vector. 
		/*  The swap method can be defined in terms of assignment.
		This requires three assignments, each of which, for a container type, is linear
		in the container's size. In a sense, then, a.swap(b) is redundant.
		This implementation guaranties a run-time complexity that is constant rather than linear.
		\params v int_vector to swap.
		Required for the Assignable Conecpt of the STL.
		*/
		void swap(int_vector &v);

		//! Resize the int_vector in terms of elements.
		/*! \param size The size to resize the int_vector in terms of elements.
		*
		*  Required for the Sequence Concept of the STL.
		*/
		void resize(const size_type size);

		//! Resize the int_vector in terms of bits.
		/*! \param size The size to resize the int_vector in terms of bits.
		*/
		void bit_resize(const size_type size);

		//! The number of elements in the int_vector.
		/*!
		Required for the Container Concept of the STL.
		\sa max_size, bit_size, capacity
		*/
		size_type size() const;

		//! Maximum size of the int_vector.
		/*!
		Required for the Container Concept of the STL.
		\sa size, bit_size, capacity
		*/
		static size_type max_size();

		//! The number of bits in the int_vector.
		/*!
		\sa size, max_size, bit_size, capacity
		*/
		size_type bit_size() const;

		//! Returns the size of the occupied bits of the int_vector.
		/*! The capacity of a int_vector is greater or equal to the
		bit_size of the vector: capacity() >= bit_size().
		\sa size, bit_size, max_size, capacity
		*/
		size_type capacity() const;

		//! Pointer to the raw data of the int_vector
		/*!
		\returns Const pointer to the raw data of the int_vector
		*/
		const uint64_t* data() const {
			return m_data;
		}

		//! Get the integer value of the binary string of length len starting at position idx in the int_vector.
		/*! \param idx Starting index of the binary representation of the integer.
		\param len Length of the binary representation of the integer. Default value is 64.
		\returns The integer value of the binary string of length len starting at position idx.
		\sa setInt, getBit, setBit
		*/
		const value_type getInt(size_type idx, const uint8_t len = 64) const;

		//! Set the bits from position idx to idx+len-1 to the binary representation of integer x.
		/*! The bit at position idx represents the least significant bit(lsb), and the bit at
		position idx+len-1 the most significant bit (msb) of x.
		\param idx Starting index of the binary representation of x.
		\param x   The integer to store in the int_vector.
		\param len The length used to store x in the int_vector. Default value is 64.
		\sa getInt, getBit, setBit
		*/
		void setInt(size_type idx, value_type x, const uint8_t len = 64);

		//! Returns the width of the integers which are accessed via the [] operator.
		/*! \returns The width of the integers which are accessed via the [] operator.
		\sa setIntWidth
		*/
		const uint8_t getIntWidth() const;

		//! Sets the width of the integers which are accessed via the [] operator, if fixedIntWidth equals 0.
		/*! \param intWidth New width of the integers accessed via the [] operator.
		\note This method has no effect if fixedIntWidth is in the range [1..64].
		\sa getIntWidth
		*/
		void setIntWidth(uint8_t intWidth);

		//! Serializes the int_vector to a stream.
		/*!
		* \return The number of bytes written to out.
		* \sa load
		*/
		size_type serialize(std::ostream &out) const;

		//! Load the int_vector for a stream.
		void load(std::istream &in);

		//! non const version of [] operator
		/*! \param i Index the i-th integer of length getIntWidth().
		* 	\return A reference to the i-th integer of length getIntWidth().
		*
		* Required for the STL Random Access Container Concept.
		*/
		inline reference operator[] (const size_type& i);

		//! const version of [] operator
		/*! \param i Index the i-th integer of length getIntWidth().
		*  \return The value of the i-th integer of length getIntWidth().
		*
		*  Required for the STL Random Access Container Concept.
		*/
		inline const_reference operator[](const size_type& i) const;

		// ! This operator should be used for displaying a specific range of the int_vector .
		/* ! \param range This const char array should consists of an integer <b>begin_idx</b> followed by two dots (..) followed by an integer <b>end_idx</b>. E.g. "0..10" or "3..7".
		\return A string containing the bit representation from index <b>begin_idx</b> to (inclusive!) <b>end_idx</b>.
		*/
		//		std::string operator[](const char* range) const;

		//! Assignment operator for the int_vector.
		/*! \param v The vector v which shuold be assigned
		\returns A copy of v.

		Required for the Assignable Conecpt of the STL.
		*/
		int_vector& operator=(const int_vector& v);

		//! Equality operator for two int_vectors.
		/*! Two int_vectors are equal if
		- capacities and sizes are equal and
		- IntWidths are equal and
		- the bits in the range [0..bit_size()-1] are equal.

		Required for the STL Equality Comparable Concept.
		\sa operator!=
		*/
		bool operator==(const int_vector<fixedIntWidth>& v) const;

		//! Inequality operator for two int_vectors.
		/*! Two int_vectors are not equal if
		- capacities and sizes are not equal or
		- IntWidths are not equal or
		- the bits in the range [0..bit_size()-1] are not equal.


		Required for the STL Equality Comparable Concept.
		\sa operator==
		*/

		bool operator!=(const int_vector& v) const;

		//! Less operator for two int_vectors
		/*! int_vector w is less than v if
		- w[i]==v[i] for i<j and w[j]<v[j] with j in [0, min(w.size(), v.size()) )
		- or w[i]==v[i] for all i < min(w.size(), v.size()) and w.size()<v.size().

		Required for the STL LessThan Comparable Concept.
		\sa operator>
		*/
		bool operator<(const int_vector& v) const; //TODO: Test

												   //! Greater operator for two int_vectors
												   /*! int_vector w is greater than v if
												   - w[i]==v[i] for i<j and w[j]>v[j] with j in [0, min(w.size(), v.size()) )
												   - or w[i]==v[i] for all i < min(w.size(), v.size()) and w.size()>v.size().

												   Required for the STL LessThan Comparable Concept.
												   \sa operator<
												   */
		bool operator>(const int_vector& v) const;//TODO: Test
												  //! Less or equal operator
												  /*!
												  Required for the STL LessThan Comparable Concept.
												  \sa operator>=, operator==, operator<
												  */
		bool operator<=(const int_vector& v) const;//TODO: Test
												   //! Greater of equal operator
												   /*!
												   Required for the STL LessThan Comparable Concept.
												   \sa operator<=, operator==, operator>
												   */
		bool operator>=(const int_vector& v) const;//TODO: Test

												   //! Iterator that points to the first element of the int_vector.
												   /*! Required for Container Concept of the STL.
												   * 	Complexity guaranty is O(1).
												   */
		const iterator begin();

		//! Iterator that points to the element after the last element of int_vector.
		/*! Required for Container Concept of the STL.
		*  Complexity guaranty is O(1).
		*/
		const iterator end();

		//! Const iterator that points to the first element of the int_vector.
		/*! Required for Container Concept of the STL.
		*/
		const const_iterator begin() const;

		//! Const iterator that points to the element after the last element of int_vector.
		/*! Required for Container Concept of the STL.
		*/
		const const_iterator end() const;

		//TODO: rbegin()
		//TODO: rend()

	private:
		//! Set the bit at position i to value b
		/* \param i Position of the bit to set to value b.
		\param b Value to set the bit at position i.
		\sa getBit
		*/
		void setBit(size_type i, const bool b = true);

		//! Get the bit at position i.
		/* \param i Position of the bit.
		\returns The boolean value of the bit at position i.
		\sa setBit
		*/
		const bool getBit(size_type i) const;


};
	/*!
	\example bit_vector_and_stl.e.cpp
	This is an example for applying algorithms of the stl to sdsl::bit_vector.
	\example int_vector_variable.e.cpp
	This is an example for a int_vector with no fixed int width.
	\example int_vector_vs_stl_vector.e.cpp
	This example compares the performance of the stl vector with the performance of the sdsl::int_vector.
	\example bit_vector_handling.e.cpp
	Example calls the load and serialize method for the sdsl::bit_vector.
	*/

	//! A proxy class that acts as a reference to an integer of length \p len bits in a int_vector.
template<class int_vector>
class int_vector_reference {
	private:
		typename int_vector::value_type * const m_word;
		const uint8_t m_offset;
		const uint8_t m_len; //!< Length of the integer referred to in bits.
							 //				uint64_t m_val; //!< Value of the integer referred to.
	public:
		//! Constructor for the reference class
		/*! \param word Pointer to the corresponding 64bit word in the int_vector.
		\param offset Offset to the starting bit (offset in [0..63])
		\param len length of the integer, should be v->getIntWidth()!!!
		*/
		int_vector_reference(typename int_vector::value_type *word, uint8_t offset, uint8_t len) :
			m_word(word), m_offset(offset), m_len(len) {};

		//! Assignment operator for the proxy class
		/*!
		The integer x is assign to the referenced
		position in the BitVector with the specified fixedIntWidth
		of the BitVector
		\param x 64bit integer to assign
		\return A const_reference to the assigned reference
		*/
		int_vector_reference& operator=(typename int_vector::value_type x) {
			bit_magic::writeInt(m_word, x, m_offset, m_len);
			return *this;
		};

		int_vector_reference& operator=(const int_vector_reference& x) {
			return *this = typename int_vector::value_type(x);
		};

		//! Cast the reference to a int_vector<>::value_type
		operator typename int_vector::value_type()const {
			return bit_magic::readInt(m_word, m_offset, m_len);
		}

		bool operator==(const int_vector_reference& x)const {
			return typename int_vector::value_type(*this) == typename int_vector::value_type(x);
		}

		bool operator<(const int_vector_reference& x)const {
			return typename int_vector::value_type(*this) < typename int_vector::value_type(x);
		}
};



template<class int_vector>
struct int_vector_iterator_base : public std::iterator<std::random_access_iterator_tag, typename int_vector::value_type, typename int_vector::difference_type/*ptrdiff_t*/> {
		typedef /*typename int_vector::size_type*/uint64_t				size_type;
		uint8_t m_offset;
		uint8_t m_len;

		int_vector_iterator_base(uint8_t offset, uint8_t len) :m_offset(offset), m_len(len)
		{}

		int_vector_iterator_base(const int_vector *v = NULL, size_type idx = 0)/*:m_offset(idx&0x3F), m_len(v->m_intWidth)*/ {
			m_offset = idx & 0x3F;
			if (v == NULL)
				m_len = 0;
			else
				m_len = v->m_intWidth;
		}
};

template<class int_vector>
struct int_vector_iterator : public int_vector_iterator_base<int_vector> {
		typedef int_vector_reference<int_vector> 	reference;
		typedef int_vector_iterator 				iterator;
		typedef reference*							pointer;
		typedef typename int_vector::size_type		size_type;
		typedef typename int_vector::difference_type difference_type;

		using int_vector_iterator_base<int_vector>::m_offset; // make m_offset easy usable
		using int_vector_iterator_base<int_vector>::m_len;    // make m_len easy usable

		typename int_vector::value_type *m_word;

		inline int_vector_iterator(int_vector *v, size_type idx) : int_vector_iterator_base<int_vector>(v, idx) {
			m_word = v->m_data + (idx >> 6);
		}


		inline int_vector_iterator(const int_vector_iterator<int_vector> &it) : int_vector_iterator_base<int_vector>(it.m_offset, it.m_len) {
			m_word = it.m_word;
		}

		reference operator*() const {
			return reference(m_word, m_offset, m_len);
		}

		//! Prefix increment of the Iterator
		iterator& operator++() {
			m_offset += m_len;
			if (m_offset >= 64) {
				m_offset &= 0x3F;
				++m_word;
			}
			return *this;
		}

		//! Postfix increment of the Iterator
		iterator operator++(int x) {
			int_vector_iterator it = *this;
			++(*this);
			return it;
		}

		//! Prefix decrement of the Iterator
		iterator& operator--() {
			m_offset -= m_len;
			if (m_offset >= 64) {
				m_offset &= 0x3F;
				--m_word;
			}
			return *this;
		}

		//! Postfix decrement of the Iterator
		iterator operator--(int x) {
			int_vector_iterator it = *this;
			--(*this);
			return it;
		}

		iterator& operator+=(difference_type i) {
			if (i<0)
				return *this -= (-i);
			difference_type t = i*m_len;
			m_word += (t >> 6);
			m_offset = (m_offset + (t & 0x3F)) & 0x3F;
			return *this;
		}

		iterator& operator-=(difference_type i) {
			if (i<0)
				return *this += (-i);
			difference_type t = i*m_len;
			m_word -= (t >> 6);
			m_offset = (m_offset - (t & 0x3F)) & 0x3F;
			return *this/* += -i*/;	// TODO: check this!!!
		}

		iterator operator+(difference_type i) const {
			iterator it = *this;
			return it += i;
		}

		iterator operator-(difference_type i) const {
			iterator it = *this;
			return it -= i;
		}

		reference operator[](difference_type i) const {
			return *(*this + i);
		}

		bool operator==(const int_vector_iterator &it)const {
			return it.m_word == m_word && it.m_offset == m_offset;
		}

		bool operator!=(const int_vector_iterator &it)const {
			return !(*this == it);
		}

		bool operator<(const int_vector_iterator &it)const {
			if (m_word == it.m_word)
				return m_offset < it.m_offset;
			return m_word < it.m_word;
		}

		bool operator>(const int_vector_iterator &it)const {
			if (m_word == it.m_word)
				return m_offset > it.m_offset;
			return m_word > it.m_word;
		}

		bool operator>=(const int_vector_iterator &it)const {
			return !(*this < it);
		}

		bool operator<=(const int_vector_iterator &it)const {
			return !(*this > it);
		}

};

template<class int_vector>
inline typename int_vector_iterator<int_vector>::difference_type operator-(const int_vector_iterator<int_vector> &x, const int_vector_iterator<int_vector> &y) {
		return  (((x.m_word - y.m_word) << 6) + x.m_offset - y.m_offset) / x.m_len;
};

template<class int_vector>
inline int_vector_iterator<int_vector> operator+(typename int_vector_iterator<int_vector>::difference_type n, const int_vector_iterator<int_vector>& it) {
		return it + n;
}

template<class int_vector>
struct int_vector_const_iterator : public int_vector_iterator_base<int_vector> {
		//	typedef SDSBitVector::value_type		reference;
		typedef typename int_vector::value_type		const_reference;
		typedef const typename int_vector::value_type*	pointer;
		typedef int_vector_const_iterator		const_iterator;
		typedef typename int_vector::size_type		size_type;
		typedef typename int_vector::difference_type difference_type;

		using int_vector_iterator_base<int_vector>::m_offset; // make m_offset easy usable
		using int_vector_iterator_base<int_vector>::m_len;    // make m_len easy usable

		const typename int_vector::value_type * m_word;

		int_vector_const_iterator(const int_vector *v = NULL, size_type idx = 0) :int_vector_iterator_base<int_vector>(v, idx), m_word(v->data() + (idx >> 6))
		{}

		int_vector_const_iterator(const int_vector_iterator<int_vector> &it) :int_vector_iterator_base<int_vector>(it.m_offset, it.m_len), m_word(it.m_word)
		{ }

		const_reference operator*() const {
			if (m_offset + m_len <= 64) {
				return ((*m_word) >> m_offset)&bit_magic::Li1Mask[m_len];
			}
			else {
				return ((*m_word) >> m_offset) | ((*(m_word + 1) & bit_magic::Li1Mask[(m_offset + m_len) & 0x3F]) << (64 - m_offset));
				//			return (((*m_word)>>m_offset) | (*(m_word+1))<<(64-m_offset)) & bit_magic::Li1Mask[(uint32_t)m_len];
			}
		}

		//! Prefix increment of the Iterator
		const_iterator& operator++() {
			m_offset += m_len;
			if (m_offset >= 64) {
				m_offset &= 0x3F;
				++m_word;
			}
			return *this;
		}

		//! Postfix increment of the Iterator
		const_iterator operator++(int x) {
			int_vector_const_iterator it = *this;
			++(*this);
			return it;
		}

		//! Prefix decrement of the Iterator
		const_iterator& operator--() {
			m_offset -= m_len;
			if (m_offset >= 64) {
				m_offset &= 0x3F;
				--m_word;
			}
			return *this;
		}

		//! Postfix decrement of the Iterator
		const_iterator operator--(int x) {
			int_vector_const_iterator it = *this;
			--(*this);
			return it;
		}

		const_iterator& operator+=(difference_type i) {
			if (i<0)
				return *this -= (-i);
			difference_type t = i*m_len;
			m_word += (t >> 6);
			if ((m_offset += (t & 0x3F))&~0x3F) {// if new offset >= 64
				++m_word;       // add one to the word
				m_offset &= 0x3F; // offset = offset mod 64
			}
			return *this;
		}

		const_iterator& operator-=(difference_type i) {
			if (i<0)
				return *this += (-i);
			difference_type t = i*m_len;
			m_word -= (t >> 6);
			if ((m_offset -= (t & 0x3F))&~0x3F) {// if new offset is < 0
				--m_word;
				m_offset &= 0x3F;
			}
			return *this;
		}

		const_iterator operator+(difference_type i) const {
			const_iterator it = *this;
			return it += i;
		}

		const_iterator operator-(difference_type i) const {
			const_iterator it = *this;
			return it -= i;
		}

		const_reference operator[](difference_type i) const {
			return *(*this + i);
		}

		bool operator==(const int_vector_const_iterator &it)const {
			return it.m_word == m_word && it.m_offset == m_offset;
		}

		bool operator!=(const int_vector_const_iterator &it)const {
			return !(*this == it);
		}

		bool operator<(const int_vector_const_iterator &it)const {
			if (m_word == it.m_word)
				return m_offset < it.m_offset;
			return m_word < it.m_word;
		}

		bool operator>(const int_vector_const_iterator &it)const {
			if (m_word == it.m_word)
				return m_offset > it.m_offset;
			return m_word > it.m_word;
		}

		bool operator>=(const int_vector_const_iterator &it)const {
			return !(*this < it);
		}

		bool operator<=(const int_vector_const_iterator &it)const {
			return !(*this > it);
		}

};

template<class int_vector>
inline typename int_vector_const_iterator<int_vector>::difference_type operator-(const int_vector_const_iterator<int_vector> &x, const int_vector_const_iterator<int_vector> &y) {
		return  (((x.m_word - y.m_word) << 6) + x.m_offset - y.m_offset) / x.m_len;
}

template<class int_vector>
inline int_vector_const_iterator<int_vector> operator+(typename int_vector_const_iterator<int_vector>::difference_type n, const int_vector_const_iterator<int_vector>& it) {
		return it + n;
}

	//std::ostream& operator<<(std::ostream&, const int_vector<1> &);

inline std::ostream& operator<<(std::ostream &os, const int_vector<1> &v) {
		for (int_vector<1>::const_iterator it = v.begin(), end = v.end(); it != end; ++it) {
			os << *it;
		}
		os << std::endl;
		return os;
}

inline std::ostream& operator<<(std::ostream &os, const int_vector<0> &v) {
		for (int_vector<0>::const_iterator it = v.begin(), end = v.end(); it != end; ++it) {
			os << *it;
			if (it + 1 != end) os << " ";
		}
		os << std::endl;
		return os;
}


	/*
	inline std::istream& operator>>(std::istream &in, int_vector<1> &v){
	char c;
	v.resize(0); // clear
	do{
	c = in.get();
	if(c=='0' || c=='1'){
	v.bit_resize(v.bit_size()+1);
	v.setBit(v.bit_size()-1,c=='1');
	}
	else{
	break;
	}
	}while(true);
	return in;
	}
	*/

	// ==== int_vector implemenation  ====

template<uint8_t fixedIntWidth, class size_type_class>
inline int_vector<fixedIntWidth, size_type_class>::int_vector(size_type elements, value_type default_value, uint8_t intWidth) :m_intWidth(intWidth) {
		m_data = NULL; // initialize m_data
		m_size = 0;
		int_vector_trait<fixedIntWidth, size_type_class>::setIntWidth(m_intWidth, intWidth);
		resize(elements);
		if (default_value == 0) {
			util::setZeroBits(*this);
		}
		else {
			for (iterator it = this->begin(), end = this->end(); it != end; ++it) {
				*it = default_value;
			}
		}
}

template<uint8_t fixedIntWidth, class size_type_class>
inline int_vector<fixedIntWidth, size_type_class>::int_vector(const int_vector &v) :m_intWidth(v.m_intWidth) {
		m_data = NULL; // initalize m_data
		m_size = 0;
		bit_resize(v.bit_size());
		if (v.capacity() > 0) {
			if (memcpy(m_data, v.data(), v.capacity() / 8) == NULL) {
				throw std::bad_alloc();
			}
		}
		int_vector_trait<fixedIntWidth, size_type_class>::setIntWidth(m_intWidth, v.m_intWidth);
		//	setIntWidth(v.getIntWidth());
}

template<uint8_t fixedIntWidth, class size_type_class>
int_vector<fixedIntWidth, size_type_class>& int_vector<fixedIntWidth, size_type_class>::operator=(const int_vector& v) {
		if (this != &v) {// if v is not the same object
			bit_resize(v.bit_size());
			if (v.bit_size()>0) {
				if (memcpy(m_data, v.data(), v.capacity() / 8) == NULL) {
					throw std::bad_alloc();
				}
			}
			setIntWidth(v.getIntWidth());
		}
		return *this;
}

	// Destructor
template<uint8_t fixedIntWidth, class size_type_class>
int_vector<fixedIntWidth, size_type_class>::~int_vector() {
		if (m_data != NULL)
			free(m_data); //fixed delete
}

template<uint8_t fixedIntWidth, class size_type_class>
void int_vector<fixedIntWidth, size_type_class>::swap(int_vector &v) {
		if (this != &v) {// if v and _this_ are not the same object 
			size_type	size = m_size;
			uint64_t 	*data = m_data;
			uint8_t		intWidth = m_intWidth;
			m_size = v.m_size;
			m_data = v.m_data;
			int_vector_trait<fixedIntWidth, size_type_class>::setIntWidth(m_intWidth, v.m_intWidth);
			v.m_size = size;
			v.m_data = data;
			int_vector_trait<fixedIntWidth, size_type_class>::setIntWidth(v.m_intWidth, intWidth);
		}
}

template<uint8_t fixedIntWidth, class size_type_class>
void int_vector<fixedIntWidth, size_type_class>::resize(const size_type size) {
		bit_resize(size * getIntWidth());
}


template<uint8_t fixedIntWidth, class size_type_class>
void int_vector<fixedIntWidth, size_type_class>::bit_resize(const size_type size) {
		bool do_realloc = ((size + 63) >> 6) != ((m_size + 63) >> 6);
		const size_type old_size = m_size;
		m_size = size;                       // set new size
		if (do_realloc) {
			uint64_t *data = NULL;
			//		std::cerr<<"realloc to m_size"<<m_size<<" bytes   m_data="<<m_data<<std::endl;
			data = (uint64_t*)realloc(m_data, (((m_size + 63) >> 6) << 3)); // if m_data == NULL realloc
																			// is equivalent to malloc
																			// if sizeInBytes == 0  
																			// realloc is equivalent to  
																			// free	
			assert(m_size == 0 || data != NULL);
			m_data = data;
			// initialize unreachable bits to 0
			if (m_size > old_size && bit_size() < capacity()) {//m_size>0 
				bit_magic::writeInt(m_data + (bit_size() >> 6), 0, (uint8_t)(bit_size() & 0x3F), (uint8_t)(capacity() - bit_size()));
			}
		}
}

template<uint8_t fixedIntWidth, class size_type_class>
inline void int_vector<fixedIntWidth, size_type_class>::setBit(size_type idx, const bool value) {
		if (idx >= m_size) {
			throw std::out_of_range("OUT_OF_RANGE_ERROR: int_vector::setBit(size_type idx, const bool value); idx >= size()!");
		}
		size_type int64_idx = idx >> 6;  // idx/64
		idx &= 0x3F; // idx inside a uint64_t
		if (value) {
			m_data[int64_idx] |= 1ULL << idx;
		}
		else {
			m_data[int64_idx] &= (bit_magic::All1Mask ^ (1ULL << idx));
		}
}

template<uint8_t fixedIntWidth, class size_type_class>
inline const bool int_vector<fixedIntWidth, size_type_class>::getBit(size_type idx)const {
		if (idx >= m_size) {
			throw std::out_of_range("OUT_OF_RANGE_ERROR: int_vector::getBit(size_type idx); idx >= size()!");
		}
		return m_data[idx >> 6] & (1ULL << (idx & 0x3F));
}

template<uint8_t fixedIntWidth, class size_type_class>
inline const typename int_vector<fixedIntWidth, size_type_class>::value_type int_vector<fixedIntWidth, size_type_class>::getInt(size_type idx, const uint8_t len)const {
		/*	if(idx+len > m_size){
		throw std::out_of_range("OUT_OF_RANGE_ERROR: SDSBitVector::getInt(size_type, uint8_t); idx+len > size()!");
		}
		if(len > 64){
		throw std::out_of_range("OUT_OF_RANGE_ERROR: SDSBitVector::getInt(size_type, uint8_t); len>64!");
		}
		*/
		return bit_magic::readInt(m_data + (idx >> 6), idx & 0x3F, len);
}

template<uint8_t fixedIntWidth, class size_type_class>
inline void int_vector<fixedIntWidth, size_type_class>::setInt(size_type idx, value_type x, const uint8_t len) {
		if (idx + len > m_size) {
			throw std::out_of_range("OUT_OF_RANGE_ERROR: int_vector::setInt(size_type, uint8_t); idx+len > size()!");
		}
		if (len > 64) {
			throw std::out_of_range("OUT_OF_RANGE_ERROR: int_vector::setInt(size_type, uint8_t); len>64!");
		}
		bit_magic::writeInt(m_data + (idx >> 6), x, idx & 0x3F, len);
}

template<uint8_t fixedIntWidth, class size_type_class>
typename int_vector<fixedIntWidth, size_type_class>::size_type int_vector<fixedIntWidth, size_type_class>::size() const {
		return m_size / m_intWidth;
}

template<uint8_t fixedIntWidth, class size_type_class>
typename int_vector<fixedIntWidth, size_type_class>::size_type int_vector<fixedIntWidth, size_type_class>::max_size() {
		return ((size_type)1) << (sizeof(size_type) * 8 - 6);// TODO: motivation for this expression
}

template<uint8_t fixedIntWidth, class size_type_class>
typename int_vector<fixedIntWidth, size_type_class>::size_type int_vector<fixedIntWidth, size_type_class>::bit_size() const {
		return m_size;
}

template<uint8_t fixedIntWidth, class size_type_class>
bool int_vector<fixedIntWidth, size_type_class>::empty() const {
		return m_size == 0;
}

template<uint8_t fixedIntWidth, class size_type_class>
inline typename int_vector<fixedIntWidth, size_type_class>::size_type int_vector<fixedIntWidth, size_type_class>::capacity() const {
		return ((m_size + 63) >> 6) << 6;
}

template<uint8_t fixedIntWidth, class size_type_class>
const uint8_t int_vector<fixedIntWidth, size_type_class>::getIntWidth()const {
		return m_intWidth;
}

template<uint8_t fixedIntWidth, class size_type_class>
void int_vector<fixedIntWidth, size_type_class>::setIntWidth(uint8_t width) {
		int_vector_trait<fixedIntWidth, size_type_class>::setIntWidth(m_intWidth, width); // delegate to trait function
}

template<uint8_t fixedIntWidth, class size_type_class>
inline typename int_vector<fixedIntWidth, size_type_class>::reference int_vector<fixedIntWidth, size_type_class>::operator[] (const size_type& idx) {
		size_type i = idx * m_intWidth;
		return reference(this->m_data + (i >> 6), i & 0x3F, m_intWidth);
}

	// specialized [] operator for 64 bit access.
template<>
inline int_vector<64>::reference int_vector<64>::operator[](const size_type& idx) {
		return *(this->m_data + idx);
}

	// specialized [] operator for 32 bit access.
template<>
inline int_vector<32>::reference int_vector<32>::operator[](const size_type& idx) {
		return *(((uint32_t*)(this->m_data)) + idx);
}

	// specialized [] operator for 16 bit access.
template<>
inline int_vector<16>::reference int_vector<16>::operator[](const size_type& idx) {
		return *(((uint16_t*)(this->m_data)) + idx);
}

	// specialized [] operator for 8 bit access.
template<>
inline int_vector<8>::reference int_vector<8>::operator[](const size_type& idx) {
		return *(((uint8_t*)(this->m_data)) + idx);
}

template<uint8_t fixedIntWidth, class size_type_class>
inline typename int_vector<fixedIntWidth, size_type_class>::const_reference int_vector<fixedIntWidth, size_type_class>::operator[](const size_type& idx)const {
		return getInt(idx * fixedIntWidth, fixedIntWidth);
}

template<>
inline int_vector<0>::const_reference int_vector<0>::operator[](const size_type& idx)const {
		return getInt(idx * m_intWidth, m_intWidth);
}

template<>
inline int_vector<64>::const_reference int_vector<64>::operator[](const size_type& idx)const {
		return *(this->m_data + idx);
}

template<>
inline int_vector<32>::const_reference int_vector<32>::operator[](const size_type& idx)const {
		return *(((uint32_t*)this->m_data) + idx);
}

template<>
inline int_vector<16>::const_reference int_vector<16>::operator[](const size_type& idx)const {
		return *(((uint16_t*)this->m_data) + idx);
}

template<>
inline int_vector<8>::const_reference int_vector<8>::operator[](const size_type& idx)const {
		return *(((uint8_t*)this->m_data) + idx);
}

template<>
inline int_vector<1>::const_reference int_vector<1>::operator[](const size_type& idx)const {
		return ((*(m_data + (idx >> 6))) >> (idx & 0x3F)) & 1;
}

template<uint8_t fixedIntWidth, class size_type_class>
inline const typename int_vector<fixedIntWidth, size_type_class>::iterator int_vector<fixedIntWidth, size_type_class>::begin() {
		//	return iterator(this, 0);
		return int_vector_trait<fixedIntWidth, size_type_class>::begin(this, m_data);
}

template<uint8_t fixedIntWidth, class size_type_class>
inline const typename int_vector<fixedIntWidth, size_type_class>::iterator int_vector<fixedIntWidth, size_type_class>::end() {
		//	return iterator(this, m_intWidth*(m_size/m_intWidth));
		return int_vector_trait<fixedIntWidth, size_type_class>::end(this, m_data, (m_size / m_intWidth));
}

template<uint8_t fixedIntWidth, class size_type_class>
inline const typename int_vector<fixedIntWidth, size_type_class>::const_iterator int_vector<fixedIntWidth, size_type_class>::begin()const {
		//	return const_iterator(this, 0);
		return int_vector_trait<fixedIntWidth, size_type_class>::begin(this, m_data);
}

template<uint8_t fixedIntWidth, class size_type_class>
inline const typename int_vector<fixedIntWidth, size_type_class>::const_iterator int_vector<fixedIntWidth, size_type_class>::end()const {
		//	return const_iterator(this, m_intWidth*(m_size/m_intWidth));
		return int_vector_trait<fixedIntWidth, size_type_class>::end(this, m_data, (m_size / m_intWidth));
}

template<uint8_t fixedIntWidth, class size_type_class>
bool int_vector<fixedIntWidth, size_type_class>::operator==(const int_vector<fixedIntWidth>& v)const {
		if (capacity() != v.capacity())
			return false;
		if (bit_size() != v.bit_size())
			return false;
		if (empty())
			return true;
		const uint64_t *data1 = v.data();
		const uint64_t *data2 = data();
		for (size_type i = 0; i < (capacity() >> 6) - 1; ++i) {
			if (*(data1++) != *(data2++))
				return false;
		}
		int8_t l = 64 - (capacity() - bit_size());
		return ((*data1)&bit_magic::Li1Mask[l]) == ((*data2)&bit_magic::Li1Mask[l]);
}

template<uint8_t fixedIntWidth, class size_type_class>
bool int_vector<fixedIntWidth, size_type_class>::operator<(const int_vector& v)const {
		size_type min_size = size();
		if (min_size > v.size())
			min_size = v.size();
		for (typename int_vector<fixedIntWidth, size_type_class>::const_iterator it = begin(), end = begin() + min_size, it_v = v.begin(); it != end; ++it, ++it_v) {
			if (*it == *it_v)
				continue;
			else
				return *it < *it_v;
		}
		return  size() < v.size();
}

template<uint8_t fixedIntWidth, class size_type_class>
bool int_vector<fixedIntWidth, size_type_class>::operator>(const int_vector& v)const {
		size_type min_size = size();
		if (min_size > v.size())
			min_size = v.size();
		for (typename int_vector<fixedIntWidth, size_type_class>::const_iterator it = begin(), end = begin() + min_size, it_v = v.begin(); it != end; ++it, ++it_v) {
			if (*it == *it_v)
				continue;
			else
				return *it > *it_v;
		}
		return  size() > v.size();
}

template<uint8_t fixedIntWidth, class size_type_class>
bool int_vector<fixedIntWidth, size_type_class>::operator<=(const int_vector& v)const {
		return *this == v or *this<v;
}

template<uint8_t fixedIntWidth, class size_type_class>
bool int_vector<fixedIntWidth, size_type_class>::operator>=(const int_vector& v)const {
		return *this == v or *this>v;
}

template<uint8_t fixedIntWidth, class size_type_class>
bool int_vector<fixedIntWidth, size_type_class>::operator!=(const int_vector& v)const {
		return !(*this == v);
}

template<uint8_t fixedIntWidth, class size_type_class>
typename int_vector<fixedIntWidth, size_type_class>::size_type int_vector<fixedIntWidth, size_type_class>::serialize(std::ostream &out) const {
		size_type written_bytes = 0;
		out.write((char *)&m_size, sizeof(m_size));
		written_bytes += sizeof(m_size);
		if (fixedIntWidth == 0) {
			out.write((char *)&m_intWidth, sizeof(m_intWidth));
			written_bytes += sizeof(m_intWidth);
		}
		uint64_t *p = m_data;
		const static size_type BLOCK_SIZE = (1 << 20);
		size_type idx = 0;
		while (idx + BLOCK_SIZE < (capacity() >> 6)) {
			out.write((char*)p, BLOCK_SIZE * sizeof(uint64_t));
			written_bytes += BLOCK_SIZE * sizeof(uint64_t);
			p += BLOCK_SIZE;
			idx += BLOCK_SIZE;
		}
		out.write((char*)p, ((capacity() >> 6) - idx) * sizeof(uint64_t));
		written_bytes += ((capacity() >> 6) - idx) * sizeof(uint64_t);
		return written_bytes;
}

template<uint8_t fixedIntWidth, class size_type_class>
void int_vector<fixedIntWidth, size_type_class>::load(std::istream &in) {
		size_type size, intWidth;
		in.read((char *)&size, sizeof(m_size));
		if (fixedIntWidth == 0) {
			in.read((char *)&intWidth, sizeof(m_intWidth));
		}
		else
			intWidth = fixedIntWidth;
		bit_resize(size);
		setIntWidth((uint8_t)intWidth);
		uint64_t *p = m_data;
		const static size_type BLOCK_SIZE = (1 << 20);
		size_type idx = 0;
		while (idx + BLOCK_SIZE < (capacity() >> 6)) {
			in.read((char*)p, BLOCK_SIZE * sizeof(uint64_t));
			p += BLOCK_SIZE;
			idx += BLOCK_SIZE;
		}
		in.read((char*)p, ((capacity() >> 6) - idx) * sizeof(uint64_t));
}

#endif // end file 

