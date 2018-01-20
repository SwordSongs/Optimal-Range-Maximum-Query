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
/*! \file enc_vector_prac.hpp
\brief enc_vector_prac.hpp contains the sdsl::enc_vector_prac class.
\author Simon Gog
*/
#include "vectors.hpp"
#include "elias_delta_coder.hpp"
#include "iterators.hpp"

#ifndef _ENC_VECTOR_PRAC_HPP_
#define _ENC_VECTOR_PRAC_HPP_

	template<uint8_t fixedIntWidth>
	struct enc_vector_prac_trait {
		typedef int_vector<0> int_vector_type;
	};

	template<>
	struct enc_vector_prac_trait<32> {
		typedef int_vector<32> int_vector_type;
	};

	//! A generic immutable space-saving vector class for unsigned positiv integers. It encodes each integer with its self-delimiting code and still provides constant time access.
	/*! The values of a enc_vector_prac are immutable after the constructor call. The class
	could be parametrized with a self-delimiting codes Coder and the sample denisty.
	*/
	template<class Coder = coder::elias_delta,
		uint32_t SampleDens = 8, uint8_t fixedIntWidth = 0>
		class enc_vector_prac {
		public:
			typedef uint64_t 							value_type;  	// STL Container requirement
			typedef random_access_const_iterator<enc_vector_prac> iterator;// STL Container requirement
			typedef iterator							const_iterator; // STL Container requirement
			typedef const value_type		 			reference;
			typedef const value_type 					const_reference;
			typedef const value_type*					const_pointer;
			typedef ptrdiff_t 							difference_type;// STL Container requirement
			typedef int_vector<>::size_type				size_type;		// STL Container requirement
			typedef Coder								coder;
			static  const uint32_t 						sample_dens = SampleDens;

			int_vector<0> 	m_z; 		// compressed bit stream
		private:
			typename enc_vector_prac_trait<fixedIntWidth>::int_vector_type   m_sample_vals_and_pointer;
			size_type		m_elements;    // number of elements
			uint32_t		m_sample_dens;

			// workaround function for the constructor
			void construct() {
				m_elements = 0;
				m_sample_dens = 8;
			}
			void copy(const enc_vector_prac &v);
		public:
			//! Default Constuctor
			enc_vector_prac() {
				construct();
			};
			//! Copy constructor
			/*! \param v The enc_vector_prac to copy.
			Required for the Assignable Concept of the STL
			*/
			enc_vector_prac(const enc_vector_prac &v);

			//! Constructor for a Container of positiv integers.
			/*! \param c A container of positive integers.
			\par The container is used to build the EncVector of the
			integer sequence.
			*/
			template<class Container>
			enc_vector_prac(const Container &c) {
				construct();
				init(c);
			};

			template<class Container>
			void init(const Container &c);

			//! Default Destructor
			~enc_vector_prac() {
			};

			//! The number of elements in the enc_vector_prac.
			/*!
			Required for the Container Concept of the STL.
			\sa max_size
			*/
			size_type size()const;

			//! Return the largest size that this container can ever have.
			/*! Required for the Container Concept of the STL.
			*/
			static size_type max_size();

			//!	Returns if the enc_vector_prac is empty.
			/*! Equivalent to size() == 0.
			*
			* 	Required for the STL Container Concept.
			*  \sa size()
			*/
			bool empty() const;

			//! Swap method for enc_vector_prac
			/*! The swap method can be defined in terms of assignment.
			This requires three assignments, each of which, for a container type, is linear
			in the container's size. In a sense, then, a.swap(b) is redundant.
			This implementation guaranties a run-time complexity that is constant rather than linear.
			\param v enc_vector_prac to swap.

			Required for the Assignable Conecpt of the STL.
			*/
			void swap(enc_vector_prac &v);

			//! Iterator that points to the first element of the enc_vector_prac.
			/*!
			* 	Required for the Container Concept of the STL.
			*  \sa end()
			*/
			const const_iterator begin()const;

			//! Iterator that points to the position after the last element of the enc_vector_prac.
			/*!
			*	Required for the Container Concept of the STL
			*  \sa begin()
			*/
			const const_iterator end()const;

			// Iterator that points to the last element of the enc_vector_prac.
			/*
			* 	Required for the Container Concept of the STL.
			*  \sa rend()
			*/
			//		reverse_iterator rbegin()const;

			// Iterator that points to the position before the first element of the enc_vector_prac.
			/*
			*	Required for the Container Concept of the STL
			*  \sa rbegin()
			*/
			//		reverse_iterator rend()const;

			//! []-operator
			/*! \param i Index of the value. \f$ i \in [0..size()-1]\f$.
			*
			*  Required for the STL Random Access Container Concept.
			*/
			value_type operator[](size_type i)const;

			//! Assignment Operator
			/*!
			*	Required for the Assignable Concept of the STL.
			*/
			enc_vector_prac& operator=(const enc_vector_prac &v);

			//! Equality Operator
			/*! Two enc_vector_pracs are equal if all member variables are equal
			*  (including the sample density of the enc_vector_pracs).
			*  \note If the sample density is not equal you should use
			*  SDSAlgorithm::equal_container_values to compare two enc_vector_pracs.
			*
			* 	Required for the Equality Comparable Concept of the STL.
			*  \sa operator!=
			*/
			bool operator==(const enc_vector_prac &v)const;

			//! Unequality Operator
			/*! Two enc_vector_pracs are unuequal if not all member variables are equal
			*  (including the sample density of the enc_vector_pracs).
			*  \note If the sample density is not equal you should use
			*  SDSAlgorithm::equal_container_values to compare two enc_vector_pracs.
			*
			* 	Required for the Equality Comparable Concept of the STL.
			*  \sa operator==
			*/
			bool operator!=(const enc_vector_prac &v)const;

			//! Serialzes the enc_vector_prac to a stream.
			/*! \param out Outstream to write the data structure.
			\return The number of written bytes.
			*/
			size_type serialize(std::ostream &out) const;

			//! Load the enc_vector_prac from a stream.
			void load(std::istream &in);

			//! Returns the ith sample of enc_vector_prac
			/*! \param i The index of the sample. 0 <= i < size()/SampleDens
			*  \return The value of the ith sample.
			*/
			value_type sample(const size_type i) const;

			uint32_t get_sample_dens() const;
			void set_sample_dens(const uint32_t sample_dens);
	};


	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	inline uint32_t enc_vector_prac<Coder, SampleDens, fixedIntWidth>::get_sample_dens() const {
		if (SampleDens == 0)
			return m_sample_dens;
		else
			return SampleDens;
	}

	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	inline void enc_vector_prac<Coder, SampleDens, fixedIntWidth>::set_sample_dens(const uint32_t sample_dens) {
		m_sample_dens = sample_dens;
	}

	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	inline typename enc_vector_prac<Coder, SampleDens, fixedIntWidth>::value_type enc_vector_prac<Coder, SampleDens, fixedIntWidth>::operator[](const size_type i)const {
		if (i + 1 == 0 || i >= m_elements) {
			throw std::out_of_range("OUT_OF_RANGE_ERROR: enc_vector_prac::operator[](size_type); i >= size()!");
			return 0;
		}
		size_type idx = i / get_sample_dens();
		return m_sample_vals_and_pointer[idx << 1] + Coder::decode_prefix_sum(m_z.data(), m_sample_vals_and_pointer[(idx << 1) + 1], i - SampleDens*idx);
	}

	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	inline typename enc_vector_prac<Coder, SampleDens, fixedIntWidth>::value_type enc_vector_prac<Coder, SampleDens, fixedIntWidth>::sample(const size_type i)const {
		if (i + 1 == 0 || i >= m_elements / get_sample_dens()) {
			throw std::out_of_range("OUT_OF_RANGE_ERROR: enc_vector_prac::sample[](size_type); i >= size()!");
			return 0;
		}
		return m_sample_vals_and_pointer[i << 1];
	}

	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	inline enc_vector_prac<>::size_type enc_vector_prac<Coder, SampleDens, fixedIntWidth>::size()const
	{
		return m_elements;
	}

	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	inline enc_vector_prac<>::size_type enc_vector_prac<Coder, SampleDens, fixedIntWidth>::max_size()
	{
		return int_vector<>::max_size() / 2; // each element could possible occupy double space with selfdelimiting codes
	}

	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	inline bool enc_vector_prac<Coder, SampleDens, fixedIntWidth>::empty()const
	{
		return 0 == m_elements;
	}


	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	void enc_vector_prac<Coder, SampleDens, fixedIntWidth>::copy(const enc_vector_prac<Coder, SampleDens, fixedIntWidth> &v) {
		m_z = v.m_z;				// copy compressed bit stream
		m_sample_vals_and_pointer = v.m_sample_vals_and_pointer;      // copy sample values
		m_elements = v.m_elements;			// copy number of stored elements
	}

	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	enc_vector_prac<Coder, SampleDens, fixedIntWidth>::enc_vector_prac(const enc_vector_prac &v) {
		copy(v);
	}

	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	enc_vector_prac<Coder, SampleDens, fixedIntWidth>& enc_vector_prac<Coder, SampleDens, fixedIntWidth>::operator=(const enc_vector_prac<Coder, SampleDens, fixedIntWidth> &v) {
		if (this != &v) {// if v and _this_ are not the same object
			copy(v);
		}
		return *this;
	}

	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	bool enc_vector_prac<Coder, SampleDens, fixedIntWidth>::operator==(const enc_vector_prac<Coder, SampleDens, fixedIntWidth> &v)const {
		if (this == &v)
			return true;
		return	 	m_elements == v.m_elements
			and	m_z == v.m_z
			and	m_sample_vals_and_pointer == v.m_sample_vals_and_pointer;
	}

	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	bool enc_vector_prac<Coder, SampleDens, fixedIntWidth>::operator!=(const enc_vector_prac<Coder, SampleDens, fixedIntWidth> &v)const {
		return !(*this == v);
	}

	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	void enc_vector_prac<Coder, SampleDens, fixedIntWidth>::swap(enc_vector_prac<Coder, SampleDens, fixedIntWidth> &v) {
		if (this != &v) {// if v and _this_ are not the same object
			m_z.swap(v.m_z);					// swap compressed bit streams
			m_sample_vals_and_pointer.swap(v.m_sample_vals_and_pointer);
			std::swap(m_elements, v.m_elements);// swap the number of elements
		}
	}


	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	template<class Container>
	void enc_vector_prac<Coder, SampleDens, fixedIntWidth>::init(const Container &c) {
		// clear BitVectors
		m_z.resize(0);
		m_elements = 0;
		//	m_inc_start.resize(0);
		m_sample_vals_and_pointer.resize(0);
		if (c.empty()) // if c is empty there is nothing to do...
			return;
		typename Container::const_iterator	it = c.begin(), end = c.end();
		typename Container::value_type 		v1 = *it, v2, max_value = 0, max_sample_value = 0, x;
		size_type samples = 0;
		size_type z_size = 0;
		for (size_type i = 0, no_sample = 0; it != end; ++it, ++i, --no_sample) {
			v2 = *it;
			if (!no_sample) { // add a sample
				no_sample = get_sample_dens();
				if (max_sample_value < v2) max_sample_value = v2;
				++samples;
			}
			else {
				if (max_value < v2 - v1) max_value = v2 - v1;
				z_size += Coder::encoding_length(v2 - v1);
			}
			v1 = v2;
		}
		//std::cerr<<"Calculate delta"<<std::endl;
		{
			//		int_vector<> delta_c( c.size()-samples, 0, sizeof(typename Container::value_type)*8 ); // Vector for difference encoding of c
			if (max_sample_value > z_size + 1)
				m_sample_vals_and_pointer.setIntWidth(bit_magic::l1BP(max_sample_value) + 1);
			else
				m_sample_vals_and_pointer.setIntWidth(bit_magic::l1BP(z_size + 1) + 1);
			m_sample_vals_and_pointer.resize(2 * samples + 2); // add 2 for last entry
															   //		int_vector<0>::iterator d_it = delta_c.begin();
															   //		int_vector<0>::iterator sv_it = m_sample_vals_and_pointer.begin();
			typename enc_vector_prac_trait<fixedIntWidth>::int_vector_type::iterator sv_it = m_sample_vals_and_pointer.begin();
			z_size = 0;
			size_type no_sample = 0;
			for (it = c.begin(); it != end; ++it, --no_sample) {
				v2 = *it;
				if (!no_sample) { // add a sample
					no_sample = get_sample_dens();
					*sv_it = v2; ++sv_it;
					*sv_it = z_size; ++sv_it;
				}
				else {
					//				*d_it = v2-v1;
					x = v2 - v1;
					//std::cerr<<"*d_it="<<*d_it<<std::endl;
					if (v2 == v1) {
						throw std::logic_error("enc_vector_prac cannot decode adjacent equal values!");
					}
					z_size += Coder::encoding_length(x);
					//				 z_size += Coder::encoding_length( *d_it );
					//				 ++d_it;
				}
				v1 = v2;
			}
			*sv_it = 0; ++sv_it;        // initialize 
			*sv_it = z_size + 1; ++sv_it; // last entry

			m_z.bit_resize(z_size);
			uint64_t *z_data = Coder::raw_data(m_z);
			uint8_t offset = 0;
			no_sample = 0;
			for (it = c.begin(); it != end; ++it, --no_sample) {
				v2 = *it;
				if (!no_sample) { // add a sample
					no_sample = get_sample_dens();
				}
				else {
					Coder::encode(v2 - v1, z_data, offset);
				}
				v1 = v2;
			}

			//		Coder::encode(delta_c, m_z); // encode delta_c to m_z
		}
		//	delta_c.resize(0);
		//std::cerr<<"Calc rank"<<std::endl;
		//std::cerr<<"Calc select"<<std::endl;
		//std::cerr<<"Finished "<<std::endl;,
		m_elements = c.size();
	}

	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	enc_vector_prac<>::size_type enc_vector_prac<Coder, SampleDens, fixedIntWidth>::serialize(std::ostream &out) const {
		size_type written_bytes = 0;
		out.write((char *)&m_elements, sizeof(m_elements));
		written_bytes += sizeof(m_elements);
		written_bytes += m_z.serialize(out);
		written_bytes += m_sample_vals_and_pointer.serialize(out);
		return written_bytes;
	}

	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	void enc_vector_prac<Coder, SampleDens, fixedIntWidth>::load(std::istream &in) {
		in.read((char *)&m_elements, sizeof(m_elements));
		m_z.load(in);
		m_sample_vals_and_pointer.load(in);
	}

	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	const typename enc_vector_prac<Coder, SampleDens, fixedIntWidth>::const_iterator enc_vector_prac<Coder, SampleDens, fixedIntWidth>::begin()const {
		return const_iterator(this, 0);
	}

	template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	const typename enc_vector_prac<Coder, SampleDens, fixedIntWidth>::const_iterator enc_vector_prac<Coder, SampleDens, fixedIntWidth>::end()const {
		return const_iterator(this, this->m_elements);
	}

#endif
