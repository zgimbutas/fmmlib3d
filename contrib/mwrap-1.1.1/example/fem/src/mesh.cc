/*
 * mesh.cc
 *   Finite element mesh implementation.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */
#include "mesh.h"

#define ME Mesh


EType::~EType()
{
}


ME::~ME()
{
    for (vector<EType*>::iterator i = materials_owned.begin();
	 i != materials_owned.end(); ++i) {
	delete(*i);
    }
}


int ME::add_node(double* x)
{
    int id = numnp();
    for (int i = 0; i < ndm_; ++i)
	X.push_back(x[i]);
    return id;
}


int ME::add_element(EType* etype, int* nodes, int nen)
{
    int id = numelt();
    int i = 0;
    for (; i < nen && i < maxnen_; ++i)  IX.push_back(nodes[i]);
    for (; i < maxnen_;            ++i)  IX.push_back(-1);
    elements.push_back(etype);
    return id;
}


void ME::add_material(EType* material)
{
    materials_owned.push_back(material);
}


void ME::set_ur(const double* ur)
{
    for (int i = 0; i < ID.size(); ++i) {
        if (ID[i] >= 0)
            U[i] = ur[ID[i]];
        if (BC[i])
            U[i] = BV[i];
    }
}


void ME::get_ur(double* ur)
{
    for (int i = 0; i < ID.size(); ++i)
        if (ID[i] >= 0)
            ur[ID[i]] = U[i];
}


void ME::get_fr(double* fr)
{
    for (int i = 0; i < ID.size(); ++i) {
        if (ID[i] >= 0)
            fr[ID[i]] = R[i];
        if (!BC[i])
            fr[ID[i]] -= BV[i];
    }
}


int ME::initialize()
{
    ID.resize(maxndf_ * numnp());
    R.resize (maxndf_ * numnp());
    U.resize (maxndf_ * numnp());
    BC.resize(maxndf_ * numnp());
    BV.resize(maxndf_ * numnp());
    return assign_ids();
}


int ME::assign_ids()
{
    numid_ = 0;
    int n = numelt();

    for (int i = 0; i < ID.size(); ++i)
	ID[i] = 0;

    for (int i = 0; i < n; ++i)
	if (elements[i])
	    elements[i]->assign_ids(this, i);

    for (int i = 0; i < ID.size(); ++i)
	ID[i] = ( (ID[i] && !BC[i]) ? numid_++ : -1 );

    return numid_;
}


void ME::assemble_F()
{
    int n = numelt();
    for (int i = 0; i < n; ++i)
	if (elements[i])
	    elements[i]->assemble_f(this, i);
}


void ME::assemble_K(MatrixAssembler* K_assembler)
{
    int n = numelt();
    for (int i = 0; i < n; ++i)
	if (elements[i])
	    elements[i]->assemble_K(this, i, K_assembler);
}
