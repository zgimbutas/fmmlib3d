/*
 * scalar1d.h
 *   Element type for 1D Laplacian.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#ifndef SCALAR1D_H
#define SCALAR1D_H

#include "etype.h"
#include "mesh.h"

class Scalar1D : public EType {
public:
    Scalar1D(double k) : k(k) {}

    void assign_ids(Mesh* mesh, int eltid);
    void assemble_f(Mesh* mesh, int eltid);
    void assemble_K(Mesh* mesh, int eltid,
                    MatrixAssembler* K_assembler);

private:
    double k;
};

#endif /* SCALAR1D_H */
