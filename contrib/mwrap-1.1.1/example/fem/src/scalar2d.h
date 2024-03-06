/*
 * scalar2d.h
 *   Element type for 2D Laplacian.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#ifndef SCALAR2D_H
#define SCALAR2D_H

#include "etype.h"
#include "mesh.h"
#include "feshapes.h"


class Scalar2D : public EType {
public:
    Scalar2D(double k) : k(k) {}

    void assign_ids(Mesh* mesh, int eltid);
    void assemble_f(Mesh* mesh, int eltid);
    void assemble_K(Mesh* mesh, int eltid,
                    MatrixAssembler* K_assembler);

private:
    double k;
    int id[4];
    void get_quad(FEShapes& quad, Mesh* mesh, int eltid);
};

#endif /* SCALAR2D_H */
