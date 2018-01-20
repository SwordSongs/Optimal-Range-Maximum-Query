#ifndef _RMQ_OPTIMAL_HPP_
#define _RMQ_OPTIMAL_HPP_

#include "balanced_parentheses_support.hpp"
#include "vectors.hpp"
#include "bitmagic.hpp"

#include <climits>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;
//using namespace sdsl;

typedef int DT;                 // type of the input data
typedef unsigned long DTidx;    // for indexing in arrays
typedef unsigned char uchar;    // must be big enough to hold block indices

class RMQ_succinct_engineered; // forward declaration

							   /* Implements the systematic RMQ-algorithm described by Fischer (2008, unpublished).
							   */
class RMQ_optimal {
public:
	// returns RMQ[i,j]
	virtual DTidx query(DTidx, DTidx);

	// returns space occupied by the index (in bytes)
	virtual unsigned long getSize();

	// return bp.excess(E[i]+(i*s)), the value of the min in the i'th block
	virtual unsigned long A(DTidx i);

	RMQ_optimal(DT* a, DTidx n);

	~RMQ_optimal();

protected:
	// size of input array
	DTidx n;

	// the DFUDS of the 2d-Min-Heap
	bit_vector U;

	// operations on U (rank/select/excess/findopen)
	balanced_parentheses_support_simple<nearest_neighbour_dictionary<30>, rank_support_v<>, select_support_mcl<0> > bp;

	// positions of block minima (of excess-values in U)
	int_vector<> E;

	// block size in U
	DTidx s;

	// number of blocks in U (always 2n+2/s)
	DTidx nb;

	// precomputed +-1RMQs
	uchar** P;

	// precomputed mask for type-calculation (lower s bits set)
	uint64_t mask;

	// method for answering +-1RMQs, assumes i <= j
	virtual DTidx RMQ_excess(DTidx i, DTidx j);

	// normal RMQs on block minima of excess sequence:
	RMQ_succinct_engineered* R;

	// count allocated memory:
	unsigned long mem;
};


/**
* Implements a type-A algorithm for RMQ, needed to answer
* RMQs on block minima of the excess sequence of the DFUDS.
* There are blocks of size s' and superblocks of size s''.
* Queries lying completely inside the blocks are handled naively.
**/
class RMQ_succinct_engineered {
public:
	// return RMQ[i,j]
	virtual DTidx query(DTidx, DTidx);

	// returns space occupied by the index (in bytes)
	virtual unsigned long getSize();

	// prepare excess sequence of U.excess(E[i]+(i*s))
	RMQ_succinct_engineered(RMQ_optimal*, DTidx);

	~RMQ_succinct_engineered();

protected:
	// number of elements
	DTidx n;

	// pointer to calling type-p scheme
	RMQ_optimal* R;

	// table M for the out-of-block queries (contains indices of block-minima)
	uchar** M;

	// because M just stores offsets (rel. to start of block), this method
	// re-calculates the true index:
	inline DTidx m(DTidx k, DTidx block) { return M[k][block] + (block*sprime); }

	// depth of table M:
	DTidx M_depth;

	// table M' for superblock-queries (contains indices of block-minima)
	DTidx** Mprime;

	// depth of table M':
	DTidx Mprime_depth;

	// block size
	DTidx sprime;

	// superblock size
	DTidx sprimeprime;

	// number of blocks (always n/sprime)
	DTidx nb;

	// number of superblocks (always n/sprimeprime)
	DTidx nsb;

	// return block-number of entry i:
	inline DTidx block(DTidx i) { return i / sprime; }

	// return superblock-number of entry i:
	inline DTidx superblock(DTidx i) { return i / sprimeprime; }

	// count allocated memory:
	unsigned long mem;

	// true if excess sequence is very small => answer queries naively
	bool sequence_very_small;
};

#endif
