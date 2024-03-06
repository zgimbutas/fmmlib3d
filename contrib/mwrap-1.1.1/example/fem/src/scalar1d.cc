/*
 * scalar1d.cc
 *   Element type for 1D Laplacian.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include "scalar1d.h"

#define ME Scalar1D


void ME::assign_ids(Mesh* mesh, int eltid) 
{
    int n1 = mesh->ix(0,eltid);
    int n2 = mesh->ix(1,eltid);
    mesh->id(0,n1) = 1;
    mesh->id(0,n2) = 1;
}


void ME::assemble_f(Mesh* mesh, int eltid)
{
    int n1 = mesh->ix(0,eltid);
    int n2 = mesh->ix(1,eltid);
    double L = mesh->x(0,n2) - mesh->x(0,n1);
    double kappa = k/L;
    double udiff = mesh->u(0,n2) - mesh->u(0,n1);
    mesh->f(0,n1) -= kappa*udiff;
    mesh->f(0,n2) += kappa*udiff;
}


void ME::assemble_K(Mesh* mesh, int eltid,
                    MatrixAssembler* K_assembler) 
{
    int n1 = mesh->ix(0,eltid);
    int n2 = mesh->ix(1,eltid);
    double L = mesh->x(0,n2) - mesh->x(0,n1);
    double kappa = k/L;
    int i[2] = {mesh->id(0,n1), mesh->id(0,n2)};
    double K_elt[4] = { kappa, -kappa, 
                        -kappa,  kappa };
    K_assembler->add_entry(i, i, K_elt, 2, 2);
}
