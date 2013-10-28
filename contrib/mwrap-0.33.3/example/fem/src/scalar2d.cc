/*
 * scalar2d.cc
 *   Element type for 2D Laplacian (4 node quad only).
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include "etype.h"
#include "mesh.h"
#include "scalar2d.h"
#include "quad2d.h"
#include "gauss2by2.h"

#define ME Scalar2D


void ME::assign_ids(Mesh* mesh, int eltid) 
{
    for (int i = 0; i < 4; ++i)
	mesh->id(0, mesh->ix(i,eltid)) = 1;
}


void ME::assemble_f(Mesh* mesh, int eltid)
{
    Quad2d quad;
    Gauss4 xi(quad);
    get_quad(quad, mesh, eltid);

    vector<double> K(4*4);
    for (xi.start(); !xi.done(); ++xi) {

        // Compute grad u at the quadrature point
        double dudx[2] = {0, 0};
        for (int j = 0; j < 4; ++j) {
            double uj = mesh->u(0,mesh->ix(j,eltid));
            dudx[0] += xi.dN(j,0) * uj;
            dudx[1] += xi.dN(j,1) * uj;
        }

        // Contribute k * dot(grad N_i, grad u) * wt
        for (int i = 0; i < 4; ++i)
            mesh->f(0,mesh->ix(i,eltid)) += 
                k * (xi.dN(i,0) * dudx[0] + 
                     xi.dN(i,1) * dudx[1]) * xi.wt();
    }
}


void ME::assemble_K(Mesh* mesh, int eltid,
		    MatrixAssembler* K_assembler) 
{
    Quad2d quad;
    Gauss4 xi(quad);
    get_quad(quad, mesh, eltid);

    vector<double> K(4*4);
    for (xi.start(); !xi.done(); ++xi)
	for (int j = 0; j < 4; ++j)
	    for (int i = 0; i < 4; ++i)
                K[4*j+i] += 
                    k * (xi.dN(i,0) * xi.dN(j,0) +
                         xi.dN(i,1) * xi.dN(j,1)) * xi.wt();
    
    K_assembler->add_entry(id, id, &(K[0]), 4, 4);
}


void ME::get_quad(FEShapes& quad, Mesh* mesh, int eltid)
{
    for (int i = 0; i < 4; ++i) {
	int n = mesh->ix(i,eltid);
	quad.set_node(i, mesh->x(n));
	id[i] = mesh->id(0,n);
    }
}
