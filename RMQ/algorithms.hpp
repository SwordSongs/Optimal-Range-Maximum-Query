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
/*! \file algorithms.hpp
\brief algorithms.hpp contains algorithms for suffixarrays.
\author Simon Gog
*/

#ifndef _ALGORITHMS_HPP
#define _ALGORITHMS_HPP

#include "algorithms_for_balanced_parentheses.hpp"
#include "algorithms_for_compressed_suffix_trees.hpp"
//#include "vectors.hpp"
#include "vectors.hpp"
#include "enc_vector_prac.hpp"
#include <iosfwd>       // forward declaration of ostream
#include <stdint.h>     // for uint64_t uint32_t declaration
#include <stdexcept>     // for exceptions
#include <iostream>
#include <cassert>
#include <stack>
#include <utility>

	// forward declarations
	//template<class Coder, uint32_t SampleDens, uint8_t fixedIntWidth>
	//class enc_vector_prac;	

template<class EncVector, class RankSupport>
class csa_sada_theo;

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
class csa_sada_prac;

template<class Csa, class Lcp, class Bp_support>
class cst_sct;





	//! A helper class containing algorithms for succinct data structures.
	/*!
	\author Simon Gog
	*/
namespace algorithm {
		//	private:
		//		algorithm(); // This helper class can not be instantiated

		// Returns if the pair (a1,a2) is lex. less or equal than (b1,b2)
template<typename I>
static bool leq(I a1, I a2, I b1, I b2) {
		return a1 < b1 || (a1 == b1 && a2 <= b2);
}

		// Returns if the triple (a1,a2,a3) is lex. less or equal than (b1,b2,b3)
template<typename I>
static bool leq(I a1, I a2, I a3, I b1, I b2, I b3) {
			return a1 < b1 || (a1 == b1 && leq(a2, a3, b2, b3));
}

		/* radixPass stable sorts the sequence a of length n  according
		to the ordering given by r. The result is written in
		b. K is the size of the alphabet of r.
		*/
template<typename I>
static void radixPass(I *a, I *b, I *r, I n, I K);

template<typename I>
static void calcSuffixArrayDC3(I *s, I *sa, I n, I K);
		//	public:	

		// Calculate the entropy of T
		double H_0(const unsigned char *c);

		// Claculate the star entropy of T, see Manzini 2001 for details
		double H_0s(const unsigned char *c);

		double H_k(const unsigned char *c, uint8_t k);

		//! Calculate the Inverse Suffix Array from a Suffix Array SA. 
		/*!
		*    - Time requirement: \f$ O( sa.size() ) \f$ i.e. linear.
		*    - Space requirement: No additional space needed.
		*  \param sa Suffix Array.
		*  \param isa Container to store the resulting inverse Suffix Array.
		*/
template<class RandomAccessContainer1, class RandomAccessContainer2>
static void sa2isa(const RandomAccessContainer1 &sa, RandomAccessContainer2 &isa);

template<class RandomAccessContainer>
static void sa2isa(const RandomAccessContainer &sa, int_vector<> &isa);

		//! Calculate the previous smaller value (psv) array for a random access container a.
		/*!
		* \param a Container to calculate the psv array.
		* \param psv Container that contains the result after the calculation.
		* \pre The array \e a contains only non negative values and a.size() == psv.size().
		* \post \f[ psv[i] = \begin{array}{rl} a.size() &\mbox{ if} \min\{ a[j] \mid j<i \} \geq a[i] \\
		max\{j\mid a[j] < a[i] \wedge j<i\} &\mbox{ otherwise.}\end{array} \f]
		*/
template<class RandomAccessContainer1, class RandomAccessContainer2>
void calculate_psv(const RandomAccessContainer1 &a, RandomAccessContainer2 &psv);


		//! Verify the result of the method caculate_psv
		/*! \return True if the RandomAccessContainer psv is the previous smaller array of a.
		*/
template<class RandomAccessContainer1, class RandomAccessContainer2>
static bool verify_psv(const RandomAccessContainer1 &a, RandomAccessContainer2 &psv);

template<class RandomAccessContainer1, class RandomAccessContainer2>
static void calculate_psv2(const RandomAccessContainer1 &a, RandomAccessContainer2 &psv);

		//! Calculate the next smaller value (nsv) array for a random access container a.
		/*!
		* \param a Container to calculate the nsv array.
		* \param nsv Container that contains the result after the calculation.
		* \pre The array \e a contains only non negative values and a.size() == nsv.size().
		* \post \f[ nsv[i] = \begin{array}{rl} 0 &\mbox{ if} \min\{ a[j] \mid j>i \} \geq a[i] \\
		min\{j\mid a[j] < a[i] \wedge j>i\} &\mbox{ otherwise.}\end{array} \f]
		*/
template<class RandomAccessContainer1, class RandomAccessContainer2>
static void calculate_nsv(const RandomAccessContainer1 &a, RandomAccessContainer2 &nsv);

		//! Verify the result of the method of calculate_nsv
template<class RandomAccessContainer1, class RandomAccessContainer2>
static bool verify_nsv(const RandomAccessContainer1 &a, RandomAccessContainer2 &nsv);

		//! Calculates the Suffix Array for a text.
		/*!
		* \param c Text (c-string) to calculate the suffix array. The lex. order is given by the ascii-codes of the characters.
		* \param len Length of the text. *(c+len)=0 and for i<len *(c+len)!=0
		* \param sa Reference to a RandomAccessContainer which will contain the result of the calculation.
		* \pre sa.size() has to be equal to len.
		*/
template<class RandomAccessContainer>
static void calculate_sa(const unsigned char *c, size_t len, RandomAccessContainer &sa);

		//! TODO: Impelement
template<class RandomAccessContainer, class Text>
static void calculate_lcp(const RandomAccessContainer &sa, const Text &text, RandomAccessContainer &lcp);

		//! Verify that a Suffix Array sa is correct for a text c
		/*!
		* \param c Text (c-string) to check the suffix array for.
		\param len Length of the text c.
		\param sa Suffix array to check.
		\return If the suffix array is correct for the text c.
		The suffix array sa is correct for c if
		- the values of sa lie in the range [0..len)
		- all values are different
		- c+sa[i] < c+sa[i+1] for all i<len-1
		*/
template<class RandomAccessContainer>
static bool verify_sa(const unsigned char *c, size_t len, const RandomAccessContainer &sa);

		//! Verify that a SelectSupport rs is correct for a int_vector<1>.
template<class SelectSupport>
inline static bool verify_select_support(const SelectSupport &ss, const int_vector<1> &b);

		//! Verify that two Containers have the same number of elements and the all corresponding (i-th) elements (0<=i<c1.size()) of both containers are equal. 
template<class Container1, class Container2>
static bool equal_container_values(const Container1 &c1, Container2 &c2);

		//! Calculate the Inverse Suffix Array inplace.
		/*! \param sa RandomAccessContainer that contains the Suffix Array and is replaced by the Inverse Suffix Array.
		* \par Time complexity
		* 		\f$ \Order{ 2*SA.size() }\f$, i.e. linear.
		* \par Space complexity
		*		Additional sa.size() bits.
		*/
template<class RandomAccessContainer>
static void sa2isa_inplace(RandomAccessContainer &sa);

		//! Calculate the \f$\Psi\f$-function for a given Burrows and Wheeler Transformation.
		/* \param bwt Burrows and Wheeler Transformation.
		* \param len Length of the bwt.
		* \param psi Container of size len for the result.
		* \par Time complexity
		*	\f$\Order{2n}\f$ i.e. linear.
		* \par Space complexity
		*	Space of bwt (\f$\Order{n}\f$ bits) + space of uncompressed \f$\Psi\f$-function (\f$\Order{4n}\f$ bits).
		*/
template<class RandomAccessContainer>
static void bwt2psi(const unsigned char *bwt, typename RandomAccessContainer::size_type len, RandomAccessContainer &psi);

		//! Calculate the \f$\Psi\f$-function for a given suffix array.
		/*! \param sa Suffix Array to calculate the \f$\Psi\f$-function for.
		* \param psi RandomAccessContainer that will contain the resulting \f$\Psi\f$-function.
		* \par Time complexity
		*		\f$ \Order{3*SA.size() }\f$, i.e. linear.
		* \par Space complexity
		* 		 Additional \f$ sa.size() \cdot \log(RandomAccessContainer::size_type)\f$ bits
		*/
template<class RandomAccessContainer1, class RandomAccessContainer2>
static void sa2psi(const RandomAccessContainer1 &sa, RandomAccessContainer2 &psi);

template<class RandomAccessContainer>
static void sa2psi(const RandomAccessContainer &sa, int_vector<> &psi);

		//! Calculate the Longest Common Prefix Table (lcptab). 
		/*! Algorithm from Kasai et al. "Linear-Time Longest-Common-Prefix Computation in Suffix Arrays and Its Applications"
		*  \param sa Suffix Array to calculate the lcptab for.
		*  \param text Text to calculate the lcptab for.
		*  \param lcp RandomAccessContainer that will contain the resulting lcptab.
		*  \par Time complexity
		*		\f$ \Order{ SA.size() } \f$, i.e. linear.
		*/
template<class RandomAccessContainer, class Text>
static void	calculate_lcp12(const RandomAccessContainer &sa, const Text &text, RandomAccessContainer& lcp);

		//! Calculate the suffix array SA out of the \f$\Psi\f$-function and \f$ SA^{-1}[0]\f$.
		/*! \param psi A \f$\Psi-\f$ function.
		* \param isa_0 \f$ SA^{-1}[0] \f$. If SA[0]=n \f$ SA^{-1}[0]=\Psi(0) \f$.
		* \param sa A RandomAccessContainer that will contain the resulting suffix array.
		* \par Time complexity
		*       \f$\Order{psi.size()}\f$, i.e. linear.
		*/
template<class RandomAccessContainer1, class RandomAccessContainer2>
static void psi2sa(const RandomAccessContainer1& psi, const typename RandomAccessContainer1::size_type isa_0, RandomAccessContainer2& sa);

template<class RandomAccessContainer>
static void psi2sa(const RandomAccessContainer& psi, const typename RandomAccessContainer::size_type isa_0, int_vector<> &sa);
		//! Calculate the inverse suffix array SA out of the \f$\Psi\f$-function and \f$ SA^{-1}[0]\f$.
		/*! \param psi A \f$\Psi-\f$ function.
		* \param isa_0 \f$ SA^{-1}[0] \f$. If SA[0]=n \f$ SA^{-1}[0]=\Psi(0) \f$.
		* \param isa A RandomAccessContainer that will contain the resulting inverse suffix array.
		* \par Time complexity
		* 	    \f$\Order{psi.size()}\f$, i.e. linear.
		*/
template<class RandomAccessContainer>
static void psi2isa(const RandomAccessContainer& psi, const typename RandomAccessContainer::size_type isa_0, RandomAccessContainer& isa);

		//! Calculates the number of occurences of the substring \f$pattern[0..pattern\_len-1]\f$ in the text of the self index.
		/*!\param index A const reference to the self-index.
		* \param pattern The substring to search for.
		* \param pattern_len The length of the substring $pattern$.
		* \return Number of occurences of pattern in the text of the self index.
		*/
template<class SelfIndex>
static typename SelfIndex::size_type count(const SelfIndex &index, typename SelfIndex::pattern_type pattern, typename SelfIndex::size_type pattern_len);

		//! Calculates the occurences of the substring \f$pattern[0..pattern\_len-1]\f$ in the text of the self index and writes the result in an \f$occ\f$.
		/*! \param index A const reference to the self-index.
		*  \param pattern The sustring to search for.
		*  \param pattern_len The length of the substring \f$pattern\f$.
		*  \param occ RandomAccessContainer in that the occurences of the substring \f$pattern[0..pattern_len-1]\f$ in the text of the self index will be written.
		*	\return Number of occurenes of the substring \f$pattern[0..pattern\_len-1]\f$ in the text of the self index.
		*/
template<class SelfIndex, class RandomAccessContainer>
static typename SelfIndex::size_type locate(const SelfIndex &index, typename SelfIndex::pattern_type pattern, typename SelfIndex::size_type pattern_len, RandomAccessContainer &occ);

		// Specializations of count and locate

		//! Calculates the count method for a (compressed) suffix tree of type Cst.
		/*! \param cst A const reference to the (compressed) suffix tree.
		*	\param pattern The pattern we seach the number of occurences.
		*	\param pattern_len The length of the pattern.
		*	\return Number of occurences of the pattern in the text of the (compressed) suffix tree.
		*/
template<class Cst>
typename Cst::size_type count(
			const Cst &cst,
			typename Cst::pattern_type pattern,
			typename Cst::size_type pattern_len);

		//! Calculates the locate method for a (compressed) suffix tree of type Cst.
		/*! \param cst A const reference to the (compressed) suffix tree.
		*	\param pattern The pattern for which we seach the occurences for.
		*	\param pattern_len The length of the pattern.
		*  \param occ A reference to a random access container in which we store the occurences of pattern in the text of the (compressed suffix array).
		*	\return Number of occurences of the pattern in the text of the (compressed) suffix tree.
		*/
template<class Cst, class RandomAccessContainer>
typename Cst::size_type locate(
			const Cst &cst,
			typename Cst::pattern_type pattern,
			typename Cst::size_type,
			RandomAccessContainer &occ
);

		//! Calculate the concatenation of edge labels from the root to the node v of the (compressed) suffix tree of type Cst.
		/*!
		*	\param cst A const reference to the (compressed) suffix tree.
		*	\param v The node where the concatenation of the edge labels ends.
		*	\return A string representing the concatenation of edge labels form the root to the node v.
		*/
template<class Cst>
std::string extract(
			const Cst &cst,
			const typename Cst::node_type &v);


		//! Calculates the count method for a compressed suffix tree of type cst_sct.
		/*!	\sa sdsl::algorithm::count
		*/
template<class Csa, class Lcp, class Bp_support>
static typename cst_sct<Csa, Lcp, Bp_support>::size_type count(
			const cst_sct<Csa, Lcp, Bp_support> &cst,
			typename cst_sct<Csa, Lcp, Bp_support>::pattern_type pattern,
			typename cst_sct<Csa, Lcp, Bp_support>::size_type pattern_len);


template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
static typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type count(
			const csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth> &csa,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::pattern_type pattern,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type pattern_len);


template<uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth, class Coder, uint32_t EncVectorSampleDens>
static typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type count(
			const csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth> &csa,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::pattern_type pattern,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type pattern_len);

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth, class RandomAccessContainer>
static typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type locate(
			const csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth> &csa,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::pattern_type pattern,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type pattern_len,
			RandomAccessContainer &occ);


template<uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth, class Coder, uint32_t EncVectorSampleDens, class RandomAccessContainer>
static typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type locate(
			const csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth> &csa,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::pattern_type pattern,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type pattern_len,
			RandomAccessContainer &occ);


template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
static typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type backward_search(
			const csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth> &csa,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::pattern_type pattern,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type pattern_len,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type &occ_begin,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type &occ_end
		);

template<uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth, class Coder, uint32_t EncVectorSampleDens>
		static typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type backward_search(
			const csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth> &csa,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::pattern_type pattern,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type pattern_len,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type &occ_begin,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type &occ_end
);

template<class EncVector, class RankSupport>
static typename csa_sada_theo<EncVector, RankSupport>::size_type backward_search(
			const csa_sada_theo<EncVector, RankSupport> &csa,
			typename csa_sada_theo<EncVector, RankSupport>::pattern_type pattern,
			typename csa_sada_theo<EncVector, RankSupport>::size_type pattern_len,
			typename csa_sada_theo<EncVector, RankSupport>::size_type &occ_begin,
			typename csa_sada_theo<EncVector, RankSupport>::size_type &occ_end
);


		//! Reconstructs the text from position \f$begin\f$ to position \f$end-1\f$ from the self index.
		/*!\param index The self-index.
		* \param begin Starting position (inclusive) of the text to extract.
		* \param end   End position (exclusive) of the text to extract.
		*/
template<class SelfIndex>
static std::string extract(const SelfIndex &index, typename SelfIndex::size_type begin, typename SelfIndex::size_type end);

template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
static std::string extract(const csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth> &csa,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type begin,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type end
);

		/*	template<class EncVector, uint32_t SampleDens, uint8_t fixedIntWidth>
		static typename csa_sada_prac<EncVector, SampleDens, fixedIntWidth>::size_type backward_search(
		const csa_sada_prac<EncVector, SampleDens, fixedIntWidth> &csa,
		typename csa_sada_prac<EncVector, SampleDens, fixedIntWidth>::pattern_type pattern,
		typename csa_sada_prac<EncVector, SampleDens, fixedIntWidth>::size_type pattern_len,
		typename csa_sada_prac<EncVector, SampleDens, fixedIntWidth>::size_type &occ_begin,
		typename csa_sada_prac<EncVector, SampleDens, fixedIntWidth>::size_type &occ_end
		);

		template<uint32_t SampleDens, uint8_t fixedIntWidth, class Coder, uint32_t EncVectorSampleDens>
		static typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, fixedIntWidth>::size_type backward_search(
		const csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, fixedIntWidth> &csa,
		typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, fixedIntWidth>::pattern_type pattern,
		typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, fixedIntWidth>::size_type pattern_len,
		typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, fixedIntWidth>::size_type &occ_begin,
		typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, fixedIntWidth>::size_type &occ_end
		);
		*/
		//	template<uint8_t w>
		//	void inplace_radix_sort(int_vector<w> &v);	
		//	template<class RandomAccessContainer, class Text>
		//	static void calulate_lcp9(const RandomAccessContainer &sa, const Text &text, RandomAccessContainer& lcp);	

template<class RandomAccessContainer, class Text, class size_type>
size_type backward_seach(RandomAccessContainer sa, Text  bwt, Text pattern, size_type l_bound, size_type r_bound);

		//};

template<class RandomAccessContainer1, class RandomAccessContainer2>
void calculate_psv(const RandomAccessContainer1 &a, RandomAccessContainer2 &psv) {
			assert(psv.size() == a.size());
			if (a.empty())
				return;
			psv[0] = psv.size();
			assert(psv[0] == psv.size());
			std::stack<typename RandomAccessContainer1::size_type> psv_index;
			typename RandomAccessContainer1::value_type min_element = a[0];
			for (typename RandomAccessContainer1::size_type i = 0; i < a.size(); ++i) {
				if (a[i] <= min_element) {
					while (!psv_index.empty())
						psv_index.pop();
					min_element = a[i];
					psv[i] = a.size();
					psv_index.push(i);
				}
				else { // a[i] > min_element => stack will not be empty
					while (a[psv_index.top()] >= a[i])
						psv_index.pop();
					psv[i] = psv_index.top();
					psv_index.push(i);
				}
			}
	}

template<class RandomAccessContainer1, class RandomAccessContainer2>
bool verify_psv(const RandomAccessContainer1 &a, RandomAccessContainer2 &psv) {
			if (a.size() != psv.size())
				return false;
			typename RandomAccessContainer1::value_type min_element = a[0];
			for (typename RandomAccessContainer1::size_type i = 0; i<a.size(); ++i) {
				if (a[i] <= min_element) {
					min_element = a[i];
					if (psv[i] != a.size())// see definition of calculate_psv
						return false;
				}
				else {
					if (psv[i] >= i)
						return false;
					if (a[psv[i]] >= a[i])
						return false;
					for (typename RandomAccessContainer1::size_type j = psv[i] + 1; j<i; ++j)
						if (a[j]<a[i])
							return false;
				}
			}
			return true;
}


template<class RandomAccessContainer1, class RandomAccessContainer2>
void calculate_psv2(const RandomAccessContainer1 &a, RandomAccessContainer2 &psv) {
			assert(psv.size() == a.size());
			if (a.empty())
				return;
			psv[0] = psv.size();
			assert(psv[0] == psv.size());
			// TODO implementing the algorithm with use of a stack
			psv[0] = psv.size();
			typedef std::pair<typename RandomAccessContainer1::value_type, typename RandomAccessContainer1::size_type> tPII;
			std::stack<tPII> psv_stack;
			typename RandomAccessContainer1::value_type min_element = a[0], ai;
			for (typename RandomAccessContainer1::size_type i = 0; i < a.size(); ++i) {
				if ((ai = a[i]) <= min_element) {
					while (!psv_stack.empty())
						psv_stack.pop();
					min_element = ai;
					psv[i] = a.size();
					psv_stack.push(tPII(ai, i));
				}
				else { // a[i] > min_element => stack will not be empty
					while (psv_stack.top().first >= ai)
						psv_stack.pop();
					psv[i] = psv_stack.top().second;
					psv_stack.push(tPII(ai, i));
				}
			}
}

		template<class RandomAccessContainer1, class RandomAccessContainer2>
		void calculate_nsv(const RandomAccessContainer1 &a, RandomAccessContainer2 &nsv) {
			assert(nsv.size() == a.size());
			if (a.empty())
				return;
			nsv[nsv.size() - 1] = 0;
			std::stack<typename RandomAccessContainer1::size_type> nsv_index;
			typename RandomAccessContainer1::value_type min_element = a[nsv.size() - 1];
			for (typename RandomAccessContainer1::size_type i = nsv.size(); i > 0; --i) {
				if (a[i - 1] <= min_element) {
					while (!nsv_index.empty())
						nsv_index.pop();
					min_element = a[i - 1];
					nsv[i - 1] = 0;
					nsv_index.push(i - 1);
				}
				else { // a[i] > min_element => stack will not be empty
					while (a[nsv_index.top()] >= a[i - 1])
						nsv_index.pop();
					nsv[i - 1] = nsv_index.top();
					nsv_index.push(i - 1);
				}
			}
		}


		template<class RandomAccessContainer1, class RandomAccessContainer2>
		bool verify_nsv(const RandomAccessContainer1 &a, RandomAccessContainer2 &nsv) {
			if (a.size() != nsv.size())
				return false;
			typename RandomAccessContainer1::value_type min_element = a[a.size() - 1];
			for (typename RandomAccessContainer1::size_type i = a.size(); i>0; --i) {
				if (a[i - 1] <= min_element) {
					min_element = a[i - 1];
					if (nsv[i - 1] != 0)// see definition of calculate_nsv
						return false;
				}
				else {
					if (nsv[i - 1] <= i - 1)
						return false;
					if (a[nsv[i - 1]] >= a[i - 1])
						return false;
					for (typename RandomAccessContainer1::size_type j = i; j<nsv[i - 1]; ++j)
						if (a[j]<a[i - 1])
							return false;
				}
			}
			return true;
		}

		/*
		template<class SelfIndex>
		typename SelfIndex::size_type algorithm::count(const SelfIndex &index, typename SelfIndex::pattern_type pattern, typename SelfIndex::size_type pattern_len){

		}
		*/

		template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
		std::string extract(const csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth> &csa,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type begin,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type end
		) {
			if (end > csa.size())
				end = csa.size();
			if (begin >= csa.size())
				begin = 0;
			if (end<begin)
				std::logic_error("algorithm::extract end < begin or (end >= csa.size() and begin>csa.size())!");
			std::string result(end - begin, ' ');
			for (typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type i = begin, order = csa(begin); i<end; ++i, order = csa.psi[order]) {
				uint16_t c_begin = 1, c_end = 257, mid;
				while (c_begin < c_end) {
					mid = (c_begin + c_end) >> 1;
					if (csa.C[mid] <= order) {
						c_begin = mid + 1;
					}
					else {
						c_end = mid;
					}
				}
				//		int res = 0;
				//		while(res<=256 and csa.C[res])
				result[i - begin] = csa.comp2char[c_begin - 1];
			}
			return result;
		}

		template<class SelfIndex>
		std::string extract(const SelfIndex &csa,
			typename SelfIndex::size_type begin,
			typename SelfIndex::size_type end
		) {
			if (end > csa.size())
				end = csa.size();
			if (begin >= csa.size())
				begin = 0;
			if (end<begin)
				std::logic_error("algorithm::extract end < begin or (end >= csa.size() and begin>csa.size())!");
			std::string result(end - begin, ' ');
			for (typename SelfIndex::size_type i = begin, order = csa(begin); i<end; ++i, order = csa.psi[order]) {
				uint16_t c_begin = 1, c_end = 257, mid;
				while (c_begin < c_end) {
					mid = (c_begin + c_end) >> 1;
					if (csa.C[mid] <= order) {
						c_begin = mid + 1;
					}
					else {
						c_end = mid;
					}
				}
				//		int res = 0;
				//		while(res<=256 and csa.C[res])
				result[i - begin] = csa.comp2char[c_begin - 1];
			}
			return result;
		}

		template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
		typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type count(
			const csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth> &csa,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::pattern_type pattern,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type pattern_len
		) {
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type occ_begin, occ_end;
			return backward_search(csa, pattern, pattern_len, occ_begin, occ_end);
		}


		template<uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth, class Coder, uint32_t EncVectorSampleDens>
		typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type count(
			const csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth> &csa,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::pattern_type pattern,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type pattern_len
		) {
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type occ_begin, occ_end;
			return backward_search(csa, pattern, pattern_len, occ_begin, occ_end);
		}

		template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth, class RandomAccessContainer>
		typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type locate(
			const csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth> &csa,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::pattern_type pattern,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type pattern_len,
			RandomAccessContainer &occ
		) {
			typedef typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type size_type;
			size_type occ_begin, occ_end, occs;
			occs = backward_search(csa, pattern, pattern_len, occ_begin, occ_end);
			occ.resize(occs);
			typename RandomAccessContainer::iterator occ_it = occ.begin();
			for (size_type i = occ_begin; i<occ_end; ++i, ++occ_it) {
				*occ_it = csa[i];
			}
			return occs;
		}


		template<uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth, class Coder, uint32_t EncVectorSampleDens, class RandomAccessContainer>
		typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type locate(
			const csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth> &csa,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::pattern_type pattern,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type pattern_len,
			RandomAccessContainer &occ
		) {
			typedef typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type size_type;
			size_type occ_begin, occ_end, occs;
			occs = backward_search(csa, pattern, pattern_len, occ_begin, occ_end);
			occ.resize(occs);
			occ.resize(occs);
			typename RandomAccessContainer::iterator occ_it = occ.begin();
			for (size_type i = occ_begin; i<occ_end; ++i, ++occ_it) {
				*occ_it = csa[i];
			}
			return occs;
		}

		template<class Cst>
		typename Cst::size_type count(
			const Cst &cst,
			typename Cst::pattern_type pattern,
			typename Cst::size_type pattern_len) {
			if (pattern_len == 0) {
				return 0;
			}
			typedef typename Cst::size_type size_type;
			typedef typename Cst::node_type node_type;
			node_type node = cst.root();
			for (size_type i = 0, char_pos = 0; cst.depth(node) < pattern_len; ++i) {
				node_type newnode = cst.child(node, (unsigned char)pattern[cst.depth(node)], char_pos);
				if (newnode == cst.root())// root node, no match found
					return 0;
				// else the first character of the newnode matches the pattern at position depth(node)
				for (size_type j = cst.depth(node) + 1; j < cst.depth(newnode) and j < pattern_len; ++j) {
					char_pos = cst.csa.psi[char_pos];
					size_type cc = cst.csa.char2comp[pattern[j]];
					if (char_pos < cst.csa.C[cc] or char_pos >= cst.csa.C[cc + 1])
						return 0;
				}
				node = newnode;
			}
			return cst.leaves_in_the_subtree(node);
		}

		template<class Cst, class RandomAccessContainer>
		typename Cst::size_type locate(
			const Cst &cst,
			typename Cst::pattern_type pattern,
			typename Cst::size_type pattern_len,
			RandomAccessContainer &occ) {
			occ.resize(0);
			typedef typename Cst::size_type size_type;
			typedef typename Cst::node_type node_type;
			node_type node = cst.root();
			for (size_type i = 0, char_pos = 0; cst.depth(node) < pattern_len; ++i) {
				node_type newnode = cst.child(node, (unsigned char)pattern[cst.depth(node)], char_pos);
				if (newnode == cst.root())// root node, no match found
					return 0;
				// else the first character of the newnode matches the pattern at position depth(node)
				for (size_type j = cst.depth(node) + 1; j < cst.depth(newnode) and j < pattern_len; ++j) {
					char_pos = cst.csa.psi[char_pos];
					size_type cc = cst.csa.char2comp[pattern[j]];
					if (char_pos < cst.csa.C[cc] or char_pos >= cst.csa.C[cc + 1])
						return 0;
				}
				node = newnode;
			}
			size_type occs = cst.leaves_in_the_subtree(node);
			occ.resize(occs);
			size_type left = cst.leftmost_suffix_array_index_in_the_subtree(node);
			size_type right = cst.rightmost_suffix_array_index_in_the_subtree(node);
			for (size_type i = left; i <= right; ++i)
				occ[i - left] = cst.csa[i];
			return occs;
		}

		template<class Cst>
		std::string extract(
			const Cst &cst,
			const typename Cst::node_type &v) {
			typedef typename Cst::size_type size_type;
			typedef typename Cst::node_type node_type;
			size_type begin;
			// first get the suffix array entry of the leftmost leaf in the subtree rooted at v
			begin = cst.csa[cst.leftmost_suffix_array_index_in_the_subtree(v)];
			// then call the extract method on the compressed suffix array
			return extract(cst.csa, begin, begin + cst.depth(v));
		}

		template<class Csa, class Lcp, class Bp_support>
		typename cst_sct<Csa, Lcp, Bp_support>::size_type count(
			const cst_sct<Csa, Lcp, Bp_support> &cst,
			typename cst_sct<Csa, Lcp, Bp_support>::pattern_type pattern,
			typename cst_sct<Csa, Lcp, Bp_support>::size_type pattern_len) {
			if (pattern_len == 0) {
				return 0;
			}
			typedef typename cst_sct<Csa, Lcp, Bp_support>::size_type size_type;
			typedef typename cst_sct<Csa, Lcp, Bp_support>::node_type node_type;
			node_type node = cst.root();
			for (size_type i = 0, char_pos = 0; node.l < pattern_len; ++i) {
				node_type newnode = cst.child(node, (unsigned char)pattern[node.l], char_pos);
				if (newnode.l == 0)// root node, no match found
					return 0;
				// else the first character of the newnode matches the pattern at position node.l
				for (size_type j = node.l + 1; j < newnode.l and j< pattern_len; ++j) {
					char_pos = cst.csa.psi[char_pos];
					size_type cc = cst.csa.char2comp[pattern[j]];
					if (char_pos < cst.csa.C[cc] or char_pos >= cst.csa.C[cc + 1])
						return 0;
				}
				node = newnode;
			}
			return cst.leaves_in_the_subtree(node);
		}


		template<uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth, class Coder, uint32_t EncVectorSampleDens>
		typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type backward_search(
			const csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth> &csa,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::pattern_type pattern,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type pattern_len,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type &occ_begin,
			typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type &occ_end
		) {
			if (pattern_len == 0) {
				occ_begin = occ_end = 0;
				return 0;
			}
			typedef typename csa_sada_prac<enc_vector_prac<Coder, EncVectorSampleDens>, SampleDens, InvSampleDens, fixedIntWidth>::size_type size_type;
			occ_begin = 0;   		// bound is inclusive
			occ_end = csa.size(); // bound is exclusive
			while (pattern_len--) {
				//		std::cout << pattern[pattern_len] << std::endl;
				size_type c = csa.char2comp[pattern[pattern_len]];
				//		std::cerr<<" c = "<< (char)csa.comp2char[c] << " occ_begin="<<occ_begin<<" occ_end="<<occ_end<<std::endl;
				size_type c_begin = csa.C[c];  // interval begin of character c (inclusive)
				size_type c_end = csa.C[c + 1];// interval end of character c (exclusive)
											   //		std::cerr<<"c_begin = "<<c_begin<<" c_end = "<<c_end<<std::endl;
											   //		std::cerr<<"occ_begin = "<<occ_begin<<" occ_end = "<<occ_end<<std::endl;
				if (c_begin == c_end) {
					occ_begin = occ_end = c_begin;
					return 0;
				}
				// logarithmic solution
				// TODO: optimize solution by only accessing sampled values of psi
				size_type begin = c_begin, end = c_end, mid;
				while (begin != end) {
					mid = (begin + end) >> 1;  // begin <= mid < end
					if (csa.psi[mid] < occ_begin) {
						begin = mid + 1;
					}
					else {
						end = mid;
					}
				}
				occ_begin = begin;
				//		assert(csa.psi[begin] >= occ_begin);
				if (begin < c_end) {
					occ_begin = begin;
				}
				else {
					return occ_begin = occ_end = 0;
				}
				end = c_end;
				if (begin >= end) {
					return occ_begin = occ_end = 0;
				}

				while (begin != end) {
					mid = (begin + end) >> 1; // begin <= mid < end
					if (csa.psi[mid] < occ_end) {
						begin = mid + 1;
					}
					else {
						end = mid;
					}
				}
				occ_end = end;
				//		std::cerr<<"occ_begin = "<<occ_begin<<" occ_end = "<<occ_end<<std::endl;
				if (occ_end <= occ_begin) {
					return occ_begin = occ_end = 0;
				}
				//		assert(m_psi[end] >= occ_end);
				if (end > c_begin) {
					occ_end = end;
				}
				else {
					return occ_end = occ_begin = 0;
				}
			}
			return occ_end - occ_begin;

		}


		template<class EncVector, uint32_t SampleDens, uint32_t InvSampleDens, uint8_t fixedIntWidth>
		typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type backward_search(
			const csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth> &csa,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::pattern_type pattern,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type pattern_len,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type &occ_begin,
			typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type &occ_end
		)
		{
			typedef typename csa_sada_prac<EncVector, SampleDens, InvSampleDens, fixedIntWidth>::size_type size_type;
			occ_begin = 0;   		// bound is inclusive
			occ_end = csa.size(); // bound is exclusive
			while (pattern_len--) {
				size_type c = csa.char2comp[pattern[pattern_len]];
				//		std::cerr<<" c = "<< (char)csa.comp2char[c] << " occ_begin="<<occ_begin<<" occ_end="<<occ_end<<std::endl;
				size_type c_begin = csa.C[c];  // interval begin of character c (inclusive)
				size_type c_end = csa.C[c + 1];// interval end of character c (exclusive)
				if (c_begin == c_end) {
					occ_begin = occ_end = c_begin;
					return 0;
				}
				// linear solution
				/*
				for(int i=c_begin; i<c_end; ++i)// search for the first position where psi[i] is greater or equal than lower bound of the old interval
				if(m_psi[i] >= occ_begin){
				occ_begin = i;
				break;
				}
				for(int i=c_end;i>c_begin;++i){// search for the last position where psi[i] is less than the upper bound of the old interval
				if(m_psi[i-1]< occ_end){
				occ_end = i;
				break;
				}
				}
				*/
				// logarithmic solution
				size_type begin = c_begin, end = c_end, mid;
				while (begin != end) {
					mid = (begin + end) >> 1;  // begin <= mid < end
					if (csa.psi[mid] < occ_begin) {
						begin = mid + 1;
					}
					else {
						end = mid;
					}
				}
				// assert(m_psi[begin] >= occ_begin);
				if (begin < c_end) {
					occ_begin = begin;
				}
				else {
					return occ_begin = occ_end = 0;
				}
				end = c_end;
				while (begin != end) {
					mid = (begin + end) >> 1; // begin <= mid < end
					if (csa.psi[mid] < occ_end) {
						begin = mid + 1;
					}
					else {
						end = mid;
					}
				}
				// assert(m_psi[end] >= occ_end);
				if (end > c_begin) {
					occ_end = end;
				}
				else {
					return occ_end = occ_begin = 0;
				}
			}
			return occ_end - occ_begin;
		}

		template<class EncVector, class RankSupport>
		typename csa_sada_theo<EncVector, RankSupport>::size_type backward_search(
			const csa_sada_theo<EncVector, RankSupport> &csa,
			typename csa_sada_theo<EncVector, RankSupport>::pattern_type pattern,
			typename csa_sada_theo<EncVector, RankSupport>::size_type pattern_len,
			typename csa_sada_theo<EncVector, RankSupport>::size_type &occ_begin,
			typename csa_sada_theo<EncVector, RankSupport>::size_type &occ_end
		) {
			typedef typename csa_sada_theo<EncVector, RankSupport>::size_type size_type;
			occ_begin = 0;   		// bound is inclusive
			occ_end = csa.size(); // bound is exclusive
			while (pattern_len--) {
				size_type c = csa.char2comp[pattern[pattern_len]];
				//		std::cerr<<" c = "<< (char)csa.comp2char[c] << " occ_begin="<<occ_begin<<" occ_end="<<occ_end<<std::endl;
				size_type c_begin = csa.C[c];  // interval begin of character c (inclusive)
				size_type c_end = csa.C[c + 1];// interval end of character c (exclusive)
				if (c_begin == c_end) {
					occ_begin = occ_end = c_begin;
					return 0;
				}
				// logarithmic solution
				size_type begin = c_begin, end = c_end, mid;
				while (begin != end) {
					mid = (begin + end) >> 1;  // begin <= mid < end
					if (csa.psi[mid] < occ_begin) {
						begin = mid + 1;
					}
					else {
						end = mid;
					}
				}
				// assert(m_psi[begin] >= occ_begin);
				if (begin < c_end) {
					occ_begin = begin;
				}
				else {
					return occ_begin = occ_end = 0;
				}
				end = c_end;
				while (begin != end) {
					mid = (begin + end) >> 1; // begin <= mid < end
					if (csa.psi[mid] < occ_end) {
						begin = mid + 1;
					}
					else {
						end = mid;
					}
				}
				// assert(m_psi[end] >= occ_end);
				if (end > c_begin) {
					occ_end = end;
				}
				else {
					return occ_end = occ_begin = 0;
				}
			}
			return occ_end - occ_begin;
		}




		template<class RandomAccessContainer>
		void bwt2psi(const unsigned char *bwt, typename RandomAccessContainer::size_type len, RandomAccessContainer &psi) {
			if (psi.size() != len)
				psi.resize(len);
			typename RandomAccessContainer::size_type C[256] = { 0 }, index_of_dollar = 0;
			for (typename RandomAccessContainer::size_type i = 0; i<len; ++i) {
				++C[bwt[i]];
				if (bwt[i] == '\0')
					index_of_dollar = i;
			}
			//	std::cerr<<"index of . = "<<index_of_dollar<<std::endl;
			for (uint16_t i = 255; i != 0; --i)
				C[i] = C[i - 1];
			C[0] = 0;
			for (uint16_t i = 1; i<256; ++i) {
				C[i] += C[i - 1];
			}
			//	assert(C[bwt[0]]==0);
			//	psi[C[bwt[0]]] = index_of_dollar;
			//	++C[bwt[0]];
			for (typename RandomAccessContainer::size_type i = 0; i<len; ++i) {
				psi[C[bwt[i]]] = i;
				++C[bwt[i]];
			}
			/*	for(typename RandomAccessContainer::size_type i=index_of_dollar+1; i<len; ++i){
			psi[C[bwt[i]]] = i;
			++C[bwt[i]];
			}
			*/
			/*
			typename RandomAccessContainer::size_type C[256] = {0}, index_of_dollar = 0;
			for(typename RandomAccessContainer::size_type i=0; i<len+1; ++i){
			++C[bwt[i]];
			if(bwt[i]=='\0')
			index_of_dollar = i;
			}
			//	std::cerr<<"index of $ = "<<index_of_dollar<<std::endl;
			for(uint16_t i=255;i!=0;--i)
			C[i] = C[i-1];
			C[0]=0;
			for(uint16_t i=1;i<256;++i){
			C[i] += C[i-1];
			}
			psi[C[bwt[0]]-1] = index_of_dollar-1;
			++C[bwt[0]];
			for(typename RandomAccessContainer::size_type i=1; i<index_of_dollar; ++i){
			psi[C[bwt[i]]-1] = i-1;
			++C[bwt[i]];
			}
			for(typename RandomAccessContainer::size_type i=index_of_dollar+1; i<len+1; ++i){
			psi[C[bwt[i]]-1] = i-1;
			++C[bwt[i]];
			}
			*/
		}

		template<class RandomAccessContainer1, class RandomAccessContainer2>
		void sa2isa(const RandomAccessContainer1 &sa, RandomAccessContainer2 &isa) {
			isa.resize(sa.size()); // init isa
			typename RandomAccessContainer1::size_type i = 0;
			for (typename RandomAccessContainer1::const_iterator sa_it = sa.begin(), end = sa.end(); sa_it != end; ++sa_it, ++i) {
				isa[*sa_it] = i;
			}
		}

		template<class RandomAccessContainer>
		void sa2isa(const RandomAccessContainer &sa, int_vector<> &isa) {
			isa.setIntWidth(bit_magic::l1BP(sa.size()) + 1);
			isa.resize(sa.size()); // init isa
			typename RandomAccessContainer::size_type i = 0;
			for (typename RandomAccessContainer::const_iterator sa_it = sa.begin(), end = sa.end(); sa_it != end; ++sa_it, ++i) {
				isa[*sa_it] = i;
			}
		}

		template<class RandomAccessContainer1, class RandomAccessContainer2>
		void sa2psi(const RandomAccessContainer1 &sa, RandomAccessContainer2 &psi) {
			RandomAccessContainer2 isa; // temporary array for the inverse suffix array
			sa2isa(sa, isa);
			psi.resize(sa.size());
			typename RandomAccessContainer1::value_type tmp; // 
			typename RandomAccessContainer2::iterator psi_it = psi.begin();
			for (typename RandomAccessContainer1::const_iterator sa_it = sa.begin(), end = sa.end(); sa_it != end; ++sa_it, ++psi_it) {
				if ((tmp = *sa_it + 1) != sa.size())
					*psi_it = isa[tmp];
				else
					*psi_it = isa[0];
			}
		}

		template<class RandomAccessContainer>
		void sa2psi(const RandomAccessContainer &sa, int_vector<> &psi) {
			int_vector<> isa; // temporary array for the inverse suffix array
			sa2isa(sa, isa);
			psi.setIntWidth(bit_magic::l1BP(sa.size()) + 1);
			psi.resize(sa.size());
			typename RandomAccessContainer::value_type tmp; // 
			int_vector<>::iterator psi_it = psi.begin();
			for (typename RandomAccessContainer::const_iterator sa_it = sa.begin(), end = sa.end(); sa_it != end; ++sa_it, ++psi_it) {
				if ((tmp = *sa_it + 1) != sa.size())
					*psi_it = isa[tmp];
				else
					*psi_it = isa[0];
			}
		}

		template<class RandomAccessContainer, class Text>
		void calculate_lcp(const RandomAccessContainer &sa, const Text &text, RandomAccessContainer &lcp) {
			lcp = sa;
			RandomAccessContainer isa;
			sa2isa(sa, isa);

			lcp[0] = 0;
			typename RandomAccessContainer::size_type i = 0, j, k, l = 0;
			for (typename RandomAccessContainer::const_iterator isa_it = isa.begin(), end = isa.end(); isa_it != end; ++isa_it, ++i) {
				if ((j = *isa_it)) {
					k = sa[j - 1];
					while (text[k + l] == text[i + l])
						++l;
					lcp[j] = l;
					l = (l == 0) ? 0 : l - 1;
				}
			}
		}

		/*
		TODO: add implementation and definition
		template<class RandomAccessContainer>
		void algorithm::calculate_lps(){

		}
		*/

		template<class RandomAccessContainer1, class RandomAccessContainer2>
		void psi2sa(const RandomAccessContainer1& psi, const typename RandomAccessContainer1::size_type isa_0, RandomAccessContainer2 &sa) {
			sa.resize(psi.size());
			if (psi.empty())
				return;
			typename RandomAccessContainer1::value_type isa_k = isa_0;
			for (typename RandomAccessContainer1::size_type k = 0, size = psi.size(); k < size; ++k, isa_k = psi[isa_k]) {
				sa[isa_k] = k;
			}
		}

		template<class RandomAccessContainer>
		void psi2sa(const RandomAccessContainer& psi, const typename RandomAccessContainer::size_type isa_0, int_vector<> &sa) {
			sa.setIntWidth(bit_magic::l1BP(psi.size()) + 1);
			sa.resize(psi.size());
			if (psi.empty())
				return;
			typename RandomAccessContainer::value_type isa_k = isa_0;
			for (typename RandomAccessContainer::size_type k = 0, size = psi.size(); k < size; ++k, isa_k = psi[isa_k]) {
				sa[isa_k] = k;
			}
		}

		template<class RandomAccessContainer>
		void psi2isa(const RandomAccessContainer& psi, const typename RandomAccessContainer::size_type isa_0, RandomAccessContainer& isa) {
			isa = psi;
			if (psi.empty())
				return;
			typename RandomAccessContainer::value_type isa_k = isa_0;
			for (typename RandomAccessContainer::size_type k = 0, size = psi.size(); k < size; ++k, isa_k = psi[isa_k]) {
				isa[k] = isa_k;
			}
		}

		template<class RandomAccessContainer>
		void calculate_sa(const unsigned char *c, size_t len, RandomAccessContainer &sa) {
			if (len == 1) { // handle special case 
				sa = RandomAccessContainer(1, 0);
				return;
			}
			uint8_t K = 0; // size of the alphabet without zero
			uint8_t occ[256] = { 0 }; // indicator if a char occurs in c
			const unsigned char *cp = c; // char pointer to c
			for (size_t i = 0; i<len; ++i) {
				occ[*(cp++)] = 1;
			}
			//	assert(occ[0]==0);// assert that '\0' does not occure in the text
			for (uint16_t i = 0; i<256; ++i) {
				if (occ[i]) ++K;
				if (i) occ[i] += occ[i - 1];
			}
			cp = c;
			if (len + 3 <= 0xFFFFFFFFULL) {
				uint32_t *s = new uint32_t[len + 3], *suffarray = new uint32_t[len];
				for (size_t i = 0; i<len; ++i) s[i] = occ[*(cp++)];

				//		for(size_t i=0; i<len;++i) std::cerr << (uint32_t)s[i]  << std::endl;
				s[len] = s[len + 1] = s[len + 2] = 0;
				calcSuffixArrayDC3(s, suffarray, (uint32_t)len, (uint32_t)K);
				for (size_t i = 0; i<len; ++i) sa[i] = suffarray[i];
				delete[] s; delete[] suffarray;
			}
			else {
				uint64_t *s = new uint64_t[len + 3], *suffarray = new uint64_t[len];
				for (size_t i = 0; i<len; ++i) s[i] = occ[(uint8_t)*(cp++)];
				s[len] = s[len + 1] = s[len + 2] = 0;
				calcSuffixArrayDC3(s, suffarray, (uint64_t)len, (uint64_t)K);
				for (size_t i = 0; i<len; ++i) sa[i] = suffarray[i];
				delete[] s; delete[] suffarray;
			}
		}

		template<typename I>
		void radixPass(I *a, I *b, I *r, I n, I K) {
			I *c = new I[K + 1]; // counter array
			for (I i = 0; i <= K; ++i) c[i] = 0;
			for (I i = 0; i < n; ++i) c[r[a[i]]]++;
			for (I i = 0, sum = 0; i <= K; ++i) {
				I t = c[i]; c[i] = sum; sum += t;
			}
			for (I i = 0; i<n; ++i) b[c[r[a[i]]]++] = a[i];
			delete[] c;
		}

		template<typename I>
		void calcSuffixArrayDC3(I *s, I *sa, I n, I K) {
			I n0 = (n + 2) / 3, n1 = (n + 1) / 3, n2 = n / 3;// number of suffixs mod {0,1,2}==0; n0>=n1>=n2
			I n02 = n0 + n2;
			I *s12 = new I[n02 + 3];  s12[n02] = s12[n02 + 1] = s12[n02 + 2] = 0;
			I *SA12 = new I[n02 + 3]; SA12[n02] = SA12[n02 + 1] = SA12[n02 + 2] = 0;
			I* s0 = new I[n0];
			I* SA0 = new I[n0];

			I j = 0;
			for (I i = 0/*, j=0*/; i < n + (n0 - n1); ++i) if (i % 3 != 0) s12[j++] = i;
			// lsb radix sort the mod 1 and mod 2 triples
			radixPass(s12, SA12, s + 2, n02, K);
			radixPass(SA12, s12, s + 1, n02, K);
			radixPass(s12, SA12, s, n02, K);

			//  cerr<<"passed radix sorts"<<endl;
			// find lexicographic names of triples
			I name = 0, c0 = 0, c1 = 0, c2 = 0;
			if (n02 > 0) c0 = s[SA12[0]] + 1;
			for (I i = 0; i < n02; ++i) {
				if (s[SA12[i]] != c0 || s[SA12[i] + 1] != c1 || s[SA12[i] + 2] != c2) {
					name++;  c0 = s[SA12[i]];  c1 = s[SA12[i] + 1];  c2 = s[SA12[i] + 2];
				}
				if (SA12[i] % 3 == 1) { s12[SA12[i] / 3] = name; } // left half
				else { s12[SA12[i] / 3 + n0] = name; } // right half
			}
			//  cerr<<"found recursive names"<<endl;
			// recurse if names are not yet unique
			if (name < n02) {
				calcSuffixArrayDC3(s12, SA12, n02, name);
				// store unique names in s12 using the suffix array 
				for (I i = 0; i < n02; ++i) s12[SA12[i]] = i + 1;
			}
			else // generate the suffix array of s12 directly
				for (I i = 0; i < n02; ++i) SA12[s12[i] - 1] = i;

			// stably sort the mod 0 suffixes from SA12 by their first character
			for (I i = 0, j = 0; i < n02; ++i) if (SA12[i] < n0) s0[j++] = 3 * SA12[i];
			radixPass(s0, SA0, s, n0, K);

			// merge sorted SA0 suffixes and sorted SA12 suffixes
			for (I p = 0, t = n0 - n1, k = 0; k < n; ++k) {
#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
				I i = GetI(); // pos of current offset 12 suffix
				I j = SA0[p]; // pos of current offset 0  suffix
				if (SA12[t] < n0 ?
					leq(s[i], s12[SA12[t] + n0], s[j], s12[j / 3]) :
					leq(s[i], s[i + 1], s12[SA12[t] - n0 + 1], s[j], s[j + 1], s12[j / 3 + n0]))
				{ // suffix from SA12 is smaller
					sa[k] = i;  t++;
					if (t == n02) { // done --- only SA0 suffixes left
						for (k++; p < n0; p++, k++) sa[k] = SA0[p];
					}
				}
				else {
					sa[k] = j;  p++;
					if (p == n0) { // done --- only SA12 suffixes left
						for (k++; t < n02; t++, k++) sa[k] = GetI();
					}
				}
			}
			delete[] s12; delete[] SA12; delete[] SA0; delete[] s0;
		}

		template<class RandomAccessContainer>
		bool verify_sa(const unsigned char *c, size_t len, const RandomAccessContainer &sa) {
			if (sa.size() != len) // check length
				return false;
			{	// check if values are in the range [0..len) and all are different
				int_vector<> occ(len);
				util::setZeroBits(occ);
				for (typename RandomAccessContainer::const_iterator it = sa.begin(), end = sa.end(); it != end; ++it) {
					size_t value = *it;
					if (value < len and !occ[value])
						occ[value] = 1;
					else
						return false;
				}
			}
			// check if the lexicographic order is correct
			if (sa.size()<2)
				return true;
			size_t v1, v2;
			v1 = v2 = sa[(size_t)0];
			for (typename RandomAccessContainer::const_iterator it = sa.begin() + 1, end = sa.end(); it != end; ++it) {
				v1 = v2;
				v2 = *it;
				size_t i = v1, j = v2;
				bool verified = false;
				for (; i != len and j != len and !verified; ++i, ++j) {
					if (c[i] < c[j])
						verified = true;
					else if (c[i]>c[j])// lex order is wrong!
						return false;
				}
				if (!verified) {
					if (i == len)
						verified = true;
					else// j==len
						return false;
				}
			}
			return true;
		}

		template<class SelectSupport>
		bool verify_select_support(const SelectSupport &ss, const int_vector<1> &b) {
			uint64_t i = 0, j = 0;
			for (int_vector<1>::const_iterator it = b.begin(), end = b.end(); it != end; ++it, ++i) {
				if (*it) {
					++j;// it's the j-th 1 detected
					if (ss.select(j) != i) return false;
				}
			}
			return true;
		}

		template<class Container1, class Container2>
		bool equal_container_values(const Container1 &c1, Container2 &c2) {
			if (c1.size() != c2.size())
				return false;
			typename Container2::const_iterator c2_it = c2.begin();
			for (typename Container1::const_iterator c1_it = c1.begin(), c1_end = c1.end(); c1_it != c1.end(); ++c1_it, ++c2_it)
				if (*c1_it != *c2_it)
					return false;
			return true;
		}
} // end namespace algorithm

#endif
