/*
 * elastic2d.h
 *   Element type for 2D elasticity
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#ifndef ELASTIC2D_H
#define ELASTIC2D_H

#include "etype.h"
#include "mesh.h"
#include "feshapes.h"


class Elastic2D : public EType {
public:
    Elastic2D(double E, double nu, const char* type = "plane strain");

    void assign_ids(Mesh* mesh, int eltid);
    void assemble_f(Mesh* mesh, int eltid);
    void assemble_K(Mesh* mesh, int eltid,
                    MatrixAssembler* K_assembler);

private:
    double D[9];
    int id[8];
    void plane_strain(double E, double nu);
    void plane_stress(double E, double nu);
    void get_quad(FEShapes& quad, Mesh* mesh, int eltid);
};

#endif /* ELASTIC2D_H */
