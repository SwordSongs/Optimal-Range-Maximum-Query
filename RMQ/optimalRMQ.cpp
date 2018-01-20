#include "optimalRMQ.hpp"

DTidx RMQ_optimal::query(DTidx ii, DTidx jj) {
	if (ii == jj) return ii;

	// transform indices because we built the 2D-Min-Heap topsy-turvy:  
	DTidx i = n - ii - 1; DTidx j = n - jj - 1;
	if (i > j) { DTidx t = i; i = j; j = t; }   // swap if necessary

	DTidx x = bp.select(i + 2);              // caution: select starts counting at 1
	DTidx y = bp.select(j + 1);
                                                                                                                       	DTidx w = RMQ_excess(x, y);
	DTidx o = bp.find_open(w);

	//cout << "x = " << x << endl;
	//cout << "y = " << y << endl;
	//cout << "w = " << w << endl;
	//cout << "bp.rank(o) : " << bp.rank(o) << endl;
	//cout << "bp.rank(w) : " << bp.rank(w) << endl;

	if (o - bp.rank(o) == i) return n - i - 1; // re-transform answer
	return n - (w - bp.rank(w)) - 1;             // re-transform answer
 }

// version with lookup:
DTidx RMQ_optimal::RMQ_excess(DTidx i, DTidx j) {
	DTidx min, min_tmp;
	DTidx b_i = i / s;    // i's block
	DTidx b_j = j / s;    // j's block
	DTidx i_pos = i - (b_i*s);     // position of i in its block
	DTidx j_pos = j - (b_j*s);     // pos. of j in its block
	uint64_t t_i = (U.data())[b_i*s / 64];     // type of i's block
	DTidx shift = (b_i % (64 / s)) * s;        // offset to start of block in word
	t_i &= (mask << shift);              // mask bits in block
	t_i >>= shift;                       // shift block to lower s bits

	if (b_i == b_j)     // only one in-block-query => look up in P
		min = (b_i*s) + P[t_i][i_pos*(s - 1) + j_pos - ((i_pos - 1)*i_pos / 2)];
	else {
		// look up minimum in i's block:
		min = (b_i*s) + P[t_i][(i_pos + 1)*(s - 1) - ((i_pos - 1)*i_pos / 2)];

		// out-of-block query:
		if (b_j > b_i + 1) {
			DTidx bm = R->query(b_i + 1, b_j - 1); // block where out-of-block-min occurs
			min_tmp = (DTidx)(E[bm] + (bm*s));          // true position of out-of-block-min
			if (bp.excess(min_tmp) < bp.excess(min)) min = min_tmp;     
		}

		// look up minimum in j's block:
		uint64_t t_j = (U.data())[b_j*s / 64]; // type of i's block
		shift = (b_j % (64 / s)) * s;          // offset to start of block in word
		t_j &= (mask << shift);              // mask bits in block
		t_j >>= shift;                       // shift block to lower s bits
		min_tmp = (b_j*s) + P[t_j][j_pos];

		if (bp.excess(min_tmp) < bp.excess(min))     min = min_tmp;    
	}

	return min;
}

unsigned long RMQ_optimal::getSize() {
	cerr << "  size of lookup table in bytes " << mem << endl;
	cerr << "  size of array storing block minima in bytes " << util::get_size_in_bytes(E) << endl;
	cerr << "  size of RMQ on E in bytes " << R->getSize() << endl;
	cerr << "  size of DFUDS in bytes " << util::get_size_in_bytes(U) << endl;
	cerr << "  size of rank/select/findopen in bytes " << util::get_size_in_bytes(bp) << endl;
	return mem +
		R->getSize() +
		util::get_size_in_bytes(E) +
		util::get_size_in_bytes(U) +
		util::get_size_in_bytes(bp);
}

RMQ_optimal::RMQ_optimal(DT* A, DTidx n) {
	cout << "preparing array of size " << n << " for RMQs..." << endl;
	this->n = n;
	mem = 0;                               // counts allocated memory
	U = bit_vector(2 * n + 2, 0);           // DFUDS of the 2d-Min-Heap
	s = sizeof(DTidx) * 8;                  // block size for construction (full word width!)
	nb = (n - 1) / s + 1;                    // number of blocks in S
	DTidx* S = new DTidx[nb];                // stack
	for (DTidx i = 0; i<nb; i++) S[i] = 0;    // empty stack
	DTidx* M = new DTidx[nb];                  // nearest block containing a 1
	M[nb - 1] = nb;                          // stopper

										    // scan input array a from left to right (in order to return leftmost min.)
										    // and build DFUDS U:
	cout << "   building DFUDS of 2D-Max-Heap..." << endl;
	DTidx U_end = 2 * n + 1;              // current end in U
	for (DT i = n - 1; i >= 0; i--) {
		U_end--; // write implicit ')'='0' at U's current end

		DTidx j = i + 1;                           // top of stack S
		DTidx bj = (j == n) ? nb : j / s;       // j's block number

		while (bj < nb && A[n - j - 1] < A[n - i - 1]) {        //修改处1  将“A[n - j - 1] > A[n - i - 1]”修改为“A[n - j - 1] < A[n - i - 1]”
			U[U_end] = 1; U_end--;        // write '('='1' at U's current end
			S[bj] &= ~(((DTidx)1) << (((bj + 1)*s) - j - 1)); // delete bit

															  // find new top of stack S:
			if (S[bj] == 0) bj = M[bj];    // jump to next block containing a '1'
			if (bj < nb) { // update j & bj
				j = ((bj + 1)*s) - bit_magic::l1BP((uint64_t)S[bj]) - 1;
				bj = j / s;
			}
		}

		// push index i on S:
		DTidx bi = i / s;                  // i's block
		if (S[bi] == 0) M[bi] = bj;      // update M
		S[bi] |= (((DTidx)1) << (((bi + 1)*s) - i - 1)); // set bit
	}

	for (DTidx i = 0; i < U_end; i++) U[i] = 1; // prepend '('s

	//***************打印U数组****************
	//cout << "DFUDS :  ";
	//for (DTidx i = 0; i < 2 * n + 2; ++i) {
	//	if (U[i] == 0)
	//		cout << ")" << " ";
	//	else
	//		cout << "(" << " ";
	//}
	//cout << endl;
	//*****************************************

	delete[] S; // free memory
	delete[] M; // ---- " ----

				// prepare U for rank_1/select0/findopen:
	cout << "   preparing DFUDS for rank/select/findopen..." << endl;
	bp = balanced_parentheses_support_simple<nearest_neighbour_dictionary<30>, rank_support_v<>, select_support_mcl<0> >(&U);

	// calc. block minima of U's excess sequence & precompute in-block queries:
	cout << "   building +-1RMQs on excess sequence..." << endl;
	s = 8;                                     // block size = 8 always best!!!
	nb = (2 * n + 1) / s + 1;                   // number of blocks in U
	mask = (1 << s) - 1;                       // set lower s bits
	DTidx nr_types = 1 << s;                     // # different block types
	DTidx nr_queries = s*(s + 1) / 2;           // # queries per block

	P = new uchar*[nr_types];                  // P: precomputed in-block queries
 	for (DTidx i = 0; i < nr_types; i++) {
		P[i] = new uchar[nr_queries];            // allocate space 
		mem += sizeof(uchar) * nr_queries;
		P[i][0] = (uchar)s;                     // init with impossible value
	}

	E = int_vector<>(nb, 0, 3);         // 3 = log2(s)，intwidth = 3，三位表示一个int型数字，8以内：000 ... 111
	DTidx bi = 0;                       // i's block
	DTidx min = 0;                      // pos. of minimum in bi

	for (DTidx i = 1; i < 2 * n + 3; i++) {          // step through DFUDS U
		if ((i%s == 0) || (i == 2 * n + 2)) {        // reached new block
			E[bi] = min - (bi*s);                    // store position of block minimum

			// precompute in-block-queries for scanned block (if necessary):
			uint64_t type = (U.data())[bi*s / 64]; // type of block
			DTidx shift = (bi % (64 / s)) * s;     // offset to start of block in word
			type &= (mask << shift);             // mask bits in block
			type >>= shift;                      // shift block to lower s bits
			if (type >= nr_types) { cout << "ALERT!!!" << endl; exit(-1); } //DEBUG
			if (P[type][0] == s) {
				DTidx p = 0;
				for (DTidx j = bi*s; j < (bi + 1)*s; j++) {
					DTidx min2 = j;
					for (DTidx k = j; k < (bi + 1)*s; k++) {
						long em = (min2 < 2 * n + 2) ? bp.excess(min2) : LONG_MAX;     
						long ek = (k < 2 * n + 2) ? bp.excess(k) : LONG_MAX;
						if (ek < em) min2 = k;             
						P[type][p++] = (uchar)(min2 - (bi*s));
					}
				}
			}

			bi++; min = i; // increase block number & init minimum
		}
		else
			if (bp.excess(i) < bp.excess(min)) min = i; // found new minimum
	}

	E[nb - 1] = min - (bi*s); // store min of last block

							  // prepare bp.excess(E[i]+(i*s)) for (normal) RMQs:
	R = new RMQ_succinct_engineered(this, nb);

	cout << "...done.\n";
}

unsigned long RMQ_optimal::A(DTidx i) {
	return bp.excess((balanced_parentheses_support_simple<nearest_neighbour_dictionary<30>, rank_support_v<1, 1>, select_support_mcl<0, 1>>::size_type)(E[i] + (i*s)));
}

RMQ_optimal::~RMQ_optimal() {
	delete R;
	for (DTidx i = 0; i < static_cast<DTidx>(1 << s); i++) delete[] P[i];
	delete[] P;
}





/************************************************************
* The following stuff is needed to support the +-1RMQs on
* the excess sequence.
************************************************************/

unsigned long RMQ_succinct_engineered::getSize() {
	return mem;
}

// adapted such that it returns leftmost min:
DTidx RMQ_succinct_engineered::query(DTidx i, DTidx j) {
	DTidx b_i = block(i);          // i's block
	DTidx b_j = block(j);          // j's block
	DTidx min, min_i, min_j, x, y; // min: to be returned
	DTidx block_difference = b_j - b_i;

	if (sequence_very_small || block_difference < 2) {
		min = i;  // i&j not far apart => scan naively
		for (x = i + 1; x <= j; x++) if (R->A(x) < R->A(min)) min = x;    
	}
	else {
		DTidx s_bi = b_i * sprime;      // start of block i
		DTidx s_bj = b_j * sprime;      // start of block j

		min = i;
		for (x = i + 1; x < s_bi + sprime; x++) if (R->A(x) < R->A(min)) min = x;    

		DTidx k, twotothek, block_tmp;  // for index calculations in M and M'
		b_i++; // block where out-of-block-query starts
		if (s_bj - s_bi - sprime < sprimeprime) { // no superblock query
			k = bit_magic::l1BP(block_difference - 2);
			twotothek = 1 << k; // 2^k
			x = m(k, b_i); y = m(k, b_j - twotothek);
			min_i = R->A(x) <= R->A(y) ? x : y;             
		}
		else { // here we have to answer a superblock-query:
			DTidx sb_i = superblock(i); // i's superblock
			DTidx sb_j = superblock(j); // j's superblock

			block_tmp = block((sb_i + 1)*sprimeprime); // end of left out-of-block-query
			k = bit_magic::l1BP(block_tmp - b_i);
			twotothek = 1 << k; // 2^k
			x = m(k, b_i); y = m(k, block_tmp + 1 - twotothek);
			min_i = R->A(x) <= R->A(y) ? x : y;                 

			if (sb_j > sb_i + 1) { // the superblock-query:
				k = bit_magic::l1BP(sb_j - sb_i - 2);
				twotothek = 1 << k;
				x = Mprime[k][sb_i + 1]; y = Mprime[k][sb_j - twotothek];
				min_j = R->A(x) <= R->A(y) ? x : y;                         
				if (R->A(min_j) < R->A(min_i)) min_i = min_j;              
			}                

			block_tmp = block(sb_j*sprimeprime); // start of right out-of-block-query
			k = bit_magic::l1BP(b_j - block_tmp);
			twotothek = 1 << k; // 2^k
			block_tmp--; // going one block to the left doesn't harm & saves tests
			x = m(k, block_tmp); y = m(k, b_j - twotothek);
			min_j = R->A(x) <= R->A(y) ? x : y;                            
			if (R->A(min_j) < R->A(min_i)) min_i = min_j;                 
		}
		if (R->A(min_i) < R->A(min)) min = min_i;                          
		for (x = s_bj; x <= j; x++) if (R->A(x) < R->A(min)) min = x;     
	}

	return min;
}

/**
* Standard Constructor.
*/
RMQ_succinct_engineered::RMQ_succinct_engineered(RMQ_optimal* R, DTidx n) {
	this->R = R;
	this->n = n;
	sprime = 1 << 6;         // block-size
	sprimeprime = 1 << 8;	   // superblock-size
	nb = block(n - 1) + 1;       // number of blocks
	nsb = superblock(n - 1) + 1; // number of superblocks

								 // The following is necessary because we've fixed s' and s'' according to the
								 // computer's word size and NOT according to the input size. This may cause
								 // the (super-)block-size to be too big, i.e., the array too small. For such
								 // small instances it isn't advisable anyway to use this data structure,
								 // because simpler methods are faster and less space consuming.
	if (nb<sprimeprime / (2 * sprime)) {
		sequence_very_small = true;
		mem = 0;
		return; // exit constructor
	}
	else sequence_very_small = false;


	// space for out-of-block- and out-of-superblock-queries:
	M_depth = (DTidx)floor(log2(((double)sprimeprime / (double)sprime)));
	M = new uchar*[M_depth];  mem = sizeof(uchar*) * M_depth;
	M[0] = new uchar[nb];     mem += sizeof(uchar) * nb;

	Mprime_depth = (DTidx)floor(log2(nsb)) + 1;
	Mprime = new DTidx*[Mprime_depth]; mem += sizeof(DTidx*) * Mprime_depth;
	Mprime[0] = new DTidx[nsb];        mem += sizeof(DTidx) * nsb;

	// fill 0'th rows of M and Mprime:
	DTidx z = 0; // minimum in current block
	DTidx q = 0; // pos. of min in current superblock
	DTidx g = 0; // number of current superblock
	DTidx p, start, end; // start & end of current block
	for (DTidx i = 0; i < nb; i++) { // step through blocks
		start = z;              // init start
		p = start;              // init minimum
		end = start + sprime;   // end of block (not inclusive!)
		if (end > n) end = n;   // last block could be smaller than sprime!
		if (R->A(z) < R->A(q)) q = z; // update minimum in superblock   

		while (++z < end) { // step through current block:
			if (R->A(z) < R->A(p)) p = z; // update minimum in block           
			if (R->A(z) < R->A(q)) q = z; // update minimum in superblock       
		}
		M[0][i] = static_cast<uchar>(p - start);                     // store index of block-minimum (offset!)
		if (z % sprimeprime == 0 || z == n) {  // reached end of superblock?
			Mprime[0][g++] = q;                  // store index of superblock-minimum
			q = z;
		}
	}

	// fill M:
	DTidx dist = 1; // always 2^(j-1)
	for (DTidx j = 1; j < M_depth; j++) {
		M[j] = new uchar[nb]; mem += sizeof(uchar) * nb;
		for (DTidx i = 0; i < nb - dist; i++) { // be careful: loop may go too far
			M[j][i] = R->A(m(j - 1, i)) <= R->A(m(j - 1, i + dist)) ?           
				M[j - 1][i] : M[j - 1][i + dist] + static_cast<uchar>(dist*sprime); // add 'skipped' elements in a
		}
		for (DTidx i = nb - dist; i < nb; i++) M[j][i] = M[j - 1][i]; // fill overhang
		dist *= 2;
	}

	// fill M':
	dist = 1; // always 2^(j-1)
	for (DTidx j = 1; j < Mprime_depth; j++) {
		Mprime[j] = new DTidx[nsb]; mem += sizeof(DTidx)*nsb;
		for (DTidx i = 0; i < nsb - dist; i++) {
			Mprime[j][i] = R->A(Mprime[j - 1][i]) <= R->A(Mprime[j - 1][i + dist]) ?    
				Mprime[j - 1][i] : Mprime[j - 1][i + dist];
		}
		for (DTidx i = nsb - dist; i < nsb; i++) Mprime[j][i] = Mprime[j - 1][i]; // overhang
		dist *= 2;
	}
}

/**
* Destructor. Deletes allocated space.
*/
RMQ_succinct_engineered::~RMQ_succinct_engineered() {
	if (sequence_very_small) return;
	for (DTidx i = 0; i < M_depth; i++) delete[] M[i];
	delete[] M;
	for (DTidx i = 0; i < Mprime_depth; i++) delete[] Mprime[i];
	delete[] Mprime;
}
