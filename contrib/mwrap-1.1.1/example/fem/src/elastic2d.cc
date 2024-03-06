/*
 * elastic2d.cc
 *   Element type for 2D elasticity (4 node quad only).
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include "etype.h"
#include "mesh.h"
#include "elastic2d.h"
#include "quad2d.h"
#include "gauss2by2.h"
#include <string.h>

#define ME Elastic2D


ME::ME(double E, double nu, const char* type)
{
    if (strcmp(type, "plane stress") == 0)
        plane_stress(E, nu);
    else
        plane_strain(E, nu);
}


void ME::plane_strain(double E, double nu)
{
    /* <generator matexpr>
    input E, nu;
    output D(3,3);

    D = E/(1+nu)/(1-2*nu) *
        [ 1-nu,   nu, 0;
            nu, 1-nu, 0;
            0,    0, (1-2*nu)/2 ];
    */
}


void ME::plane_stress(double E, double nu)
{
    /* <generator matexpr>
    input E, nu;
    output D(3,3);

    D = E/(1-nu*nu) *
        [  1, nu, 0;
          nu,  1, 0;
           0,  0, (1-nu)/2 ];
    */
}


void ME::assign_ids(Mesh* mesh, int eltid) 
{
    for (int i = 0; i < 4; ++i) {
        int ni = mesh->ix(i,eltid);
        mesh->id(0,ni) = 1;
        mesh->id(1,ni) = 1;
    }
}


void ME::assemble_f(Mesh* mesh, int eltid)
{
    Quad2d quad;
    Gauss4 xi(quad);
    get_quad(quad, mesh, eltid);

    vector<double> K(4*4);
    for (xi.start(); !xi.done(); ++xi) {

        double eps[3] = {0, 0, 0};
        for (int j = 0; j < 4; ++j) {
            const double* u = mesh->u(mesh->ix(j,eltid));
            double Nj_x = xi.dN(j,0);
            double Nj_y = xi.dN(j,1);

            /* <generator matexpr>
            // Contribute strain at quadrature point
            
            input Nj_x, Nj_y, D(3,3);
            inout eps(3);
            input u(2);

            Bj = [Nj_x,    0;
                     0, Nj_y;
                  Nj_y, Nj_x];

            eps += Bj*u;
            */
        }

        for (int i = 0; i < 4; ++i) {
            double* f = mesh->f(mesh->ix(i,eltid));
            double Ni_x = xi.dN(i,0);
            double Ni_y = xi.dN(i,1);
            double w = xi.wt();

            /* <generator matexpr>
            // Contribute B_i'*D*B(u) * w
            
            input Ni_x, Ni_y, D(3,3), w;
            input eps(3);
            inout f(2);

            Bi = [Ni_x,    0;
                     0, Ni_y;
                  Ni_y, Ni_x];

            f += Bi'*D*eps*w;
            */
        }

    }
}


void ME::assemble_K(Mesh* mesh, int eltid,
                    MatrixAssembler* K_assembler) 
{
    Quad2d quad;
    Gauss4 xi(quad);
    get_quad(quad, mesh, eltid);

    vector<double> K(8*8);
    for (xi.start(); !xi.done(); ++xi) {
        double w = xi.wt();
        for (int j = 0; j < 4; ++j) {
            double Nj_x = xi.dN(j,0);
            double Nj_y = xi.dN(j,1);
            for (int i = 0; i < 4; ++i) {
                double Ni_x = xi.dN(i,0);
                double Ni_y = xi.dN(i,1);
                double* Knodal = &(K[16*j +2*i+ 0]);

                /* <generator matexpr>
                // B-matrix based displacement formulation
                // Isotropic plane strain constitutive tensor

                input D(3,3), w;
                input Ni_x, Ni_y, Nj_x, Nj_y;
                inout Knodal[8](2,2);

                Bi = [Ni_x,    0;
                         0, Ni_y;
                      Ni_y, Ni_x];

                Bj = [Nj_x,    0;
                         0, Nj_y;
                      Nj_y, Nj_x];

                Knodal += Bi'*D*Bj * w;
                */
            }
        }
    }
    
    K_assembler->add_entry(id, id, &(K[0]), 8, 8);
}


void ME::get_quad(FEShapes& quad, Mesh* mesh, int eltid)
{
    for (int i = 0; i < 4; ++i) {
	int n = mesh->ix(i,eltid);
	quad.set_node(i, mesh->x(n));
	id[2*i+0] = mesh->id(0,n);
	id[2*i+1] = mesh->id(1,n);
    }
}
