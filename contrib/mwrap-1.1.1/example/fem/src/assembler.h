/*
 * assembler.h
 *   Compressed sparse column matrix assembler interface.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include <vector>
using std::vector;

/*
 * The MatrixAssembler class is responsible for assembling a finite element
 * matrix.  In pure MATLAB, the assembly operation looks something like
 *
 *   for i = 1:num_elements
 *     [Ke, idx] = element_stiffness(i);
 *     Ivalid    = find(idx > 0 & idx < N);
 *     idx       = idx(Ivalid);
 *     Ke        = Ke(Ivalid,Ivalid);
 *     K(idx,idx) = K(idx,idx) + Ke;
 *   end
 *
 * The method add_entry is equivalent to most of the body of this loop --
 * it removes out-of-range indices and accumulates the rest of the element
 * contribution.
 *
 * Elements in the matrix can be stored in one of two ways.  First,
 * there are elements that go into the current compressed sparse
 * column matrix structure.  Then there are elements that correspond
 * to indices that were not in the matrix the last time we built the
 * compressed sparse column structure.  When it is time to assemble
 * the matrix, we use the compress method to take all these extra
 * elements and merge them into the compressed sparse column indexing
 * structure.  This means that after K is assembled once, any we can
 * assemble other matrices with the same structure using a little less
 * time and memory.
 */

class MatrixAssembler {
public:
    MatrixAssembler(int m, int n) : m(m), n(n), jc(n+1) {}
    MatrixAssembler(int m, int n, int coord_nnz) : 
        m(m), n(n), jc(n+1), coords(coord_nnz) {}

    void add_entry(int i, int j, double Aij);
    void add_entry(const int* i, const int* j, double* Aij,
		   int m_elt, int n_elt);

    void wipe();
    void pack_cache();
    void compress();

    int     get_m()  { return m; }
    int     get_n()  { return n; }
    int*    get_jc() { return &(jc[0]); }
    int*    get_ir() { return &(ir[0]); }
    double* get_pr() { return &(pr[0]); }

    int cache_nnz() { return coords.size(); }
    int csc_nnz()   { return pr.size();     }

private:

    struct Coord {
	Coord() {}
	Coord(int i, int j, double Aij) : i(i), j(j), Aij(Aij) {}
	bool operator<(const Coord& coord) const {
	    return ((j < coord.j) || (j == coord.j && i < coord.i));
	}
	int i, j;
	double Aij;
    };

    int m, n;

    // Things that fit in the compressed sparse representation
    vector<int> jc;
    vector<int> ir;
    vector<double> pr;

    // Cache of entries that didn't fit in the compressed sparse form
    vector<Coord> coords;
};

#endif /* ASSEMBLER_H */
