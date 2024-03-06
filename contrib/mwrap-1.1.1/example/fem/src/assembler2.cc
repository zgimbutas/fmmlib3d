/*
 * assembler.cc
 *   Compressed sparse column matrix assembler implementation.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include "assembler.h"
#include <algorithm>
#include <stdio.h>

#define ME MatrixAssembler2

using std::sort;
using std::fill;


/*
 * Add A(i,j) += Aij.  If there's space in the existing compressed
 * sparse column data structure, add it there; otherwise, stash it in
 * the coords cache vector.
 */
void ME::add_entry(int i, int j, double Aij)
{
    if (i < 0 || i >= m || j < 0 || j >= n)
	return;
    for (int ii = jc[j]; ii < jc[j+1]; ++ii) {
	if (ir[ii] == i) {
	    pr[ii] += Aij;
	    return;
	}
    }
    coords.push_back(Coord(i,j,Aij));
}


/*
 * Add an element submatrix.
 */
void ME::add_entry(const int* i, const int* j, double* Aij,
		   int m_elt, int n_elt)
{
    for (int jj = 0; jj < n_elt; ++jj)
	for (int ii = 0; ii < m_elt; ++ii)
	    add_entry(i[ii], j[jj], Aij[jj*m_elt+ii]);
}


/*
 * Wipe the input matrix.  Note that the compressed sparse column index
 * structure remains intact, so we can re-assemble without doing too much
 * work.
 */
void ME::wipe()
{
    coords.resize(0);
    fill(get_pr(), get_pr()+jc[n], 0);
}


/*
 * Sort the entries in the coords cache, and merge (by summing) contributions
 * to the same position.
 */
void ME::pack_cache()
{
    if (coords.size() > 0) {
        sort(coords.begin(), coords.end());
	int i_merge = 0;
	for (int i_coord = 1; i_coord < coords.size(); ++i_coord) {
	    if (coords[i_merge] < coords[i_coord]) 
		coords[++i_merge] = coords[i_coord];
	    else 
		coords[i_merge].Aij += coords[i_coord].Aij;
	}
	coords.resize(i_merge+1);
    }
}


/*
 * Pack the coordinate cache, then merge sort the contents of the compressed
 * sparse column data structure with the contents of the entry cache.
 */
void ME::compress()
{
    pack_cache();

    ir.resize(ir.size() + coords.size());
    pr.resize(pr.size() + coords.size());

    int i_coord = coords.size()-1;  // Next input coord to process
    int i_csc   = jc[n]-1;          // Next input CSC entry to process
    int i_merge = pr.size()-1;      // Next output CSC entry to process

    jc[n] = pr.size();
    for (int j = n; j > 0; --j) {

	// Merge column from cache with column from previous CSC
	while (i_coord >= 0 && coords[i_coord].j == j-1 &&
	       i_csc >= 0   && i_csc >= jc[j-1]) {
	    if (coords[i_coord].i > ir[i_csc]) {
		ir[i_merge] = coords[i_coord].i;
		pr[i_merge] = coords[i_coord].Aij;
		--i_coord;
	    } else {
		ir[i_merge] = ir[i_csc];
		pr[i_merge] = pr[i_csc];
		--i_csc;
	    }
	    --i_merge;
	}

	// Copy stragglers from coord list
	while (i_coord >= 0 && coords[i_coord].j == j-1) {
	    ir[i_merge] = coords[i_coord].i;
	    pr[i_merge] = coords[i_coord].Aij;
	    --i_coord;
	    --i_merge;
	}

	// Copy stragglers from CSC list
	while (i_csc >= 0 && i_csc >= jc[j-1]) {
	    ir[i_merge] = ir[i_csc];
	    pr[i_merge] = pr[i_csc];
	    --i_csc;
	    --i_merge;
	}

	// Update the column count
	jc[j-1] = i_merge+1;
    }

    coords.resize(0);
}

