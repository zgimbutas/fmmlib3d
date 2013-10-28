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
    /* <generated matexpr> */ {
    double tmp1_ = E;
    double tmp2_ = nu;
    double tmp4_ = 1.0 + tmp2_;
    double tmp5_ = tmp1_ / tmp4_;
    double tmp7_ = 2.0 * tmp2_;
    double tmp8_ = 1.0 - tmp7_;
    double tmp9_ = tmp5_ / tmp8_;
    double tmp10_ = 1.0 - tmp2_;
    double tmp12_ = tmp8_ / 2.0;
    double tmp13_ = tmp9_ * tmp10_;
    double tmp14_ = tmp9_ * tmp2_;
    double tmp15_ = tmp9_ * tmp12_;
    D[0*3+0] = tmp13_;
    D[0*3+1] = tmp14_;
    D[0*3+2] = 0;
    D[1*3+0] = tmp14_;
    D[1*3+1] = tmp13_;
    D[1*3+2] = 0;
    D[2*3+0] = 0;
    D[2*3+1] = 0;
    D[2*3+2] = tmp15_;
    } /* </generated matexpr> */
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
    /* <generated matexpr> */ {
    double tmp1_ = E;
    double tmp2_ = nu;
    double tmp4_ = tmp2_ * tmp2_;
    double tmp5_ = 1.0 - tmp4_;
    double tmp6_ = tmp1_ / tmp5_;
    double tmp8_ = 1.0 - tmp2_;
    double tmp10_ = tmp8_ / 2.0;
    double tmp11_ = tmp6_ * tmp2_;
    double tmp12_ = tmp6_ * tmp10_;
    D[0*3+0] = tmp6_;
    D[0*3+1] = tmp11_;
    D[0*3+2] = 0;
    D[1*3+0] = tmp11_;
    D[1*3+1] = tmp6_;
    D[1*3+2] = 0;
    D[2*3+0] = 0;
    D[2*3+1] = 0;
    D[2*3+2] = tmp12_;
    } /* </generated matexpr> */
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
            /* <generated matexpr> */ {
            double tmp1_ = Nj_x;
            double tmp2_ = Nj_y;
            double tmp3_ = eps[0*3+0];
            double tmp4_ = eps[0*3+1];
            double tmp5_ = eps[0*3+2];
            double tmp6_ = u[0*2+0];
            double tmp7_ = u[0*2+1];
            double tmp8_ = tmp1_ * tmp6_;
            double tmp9_ = tmp2_ * tmp7_;
            double tmp10_ = tmp2_ * tmp6_;
            double tmp11_ = tmp1_ * tmp7_;
            double tmp12_ = tmp10_ + tmp11_;
            double tmp13_ = tmp3_ + tmp8_;
            double tmp14_ = tmp4_ + tmp9_;
            double tmp15_ = tmp5_ + tmp12_;
            eps[0*3+0] = tmp13_;
            eps[0*3+1] = tmp14_;
            eps[0*3+2] = tmp15_;
            } /* </generated matexpr> */
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
            /* <generated matexpr> */ {
            double tmp1_ = Ni_x;
            double tmp2_ = Ni_y;
            double tmp3_ = D[0*3+0];
            double tmp4_ = D[0*3+1];
            double tmp5_ = D[0*3+2];
            double tmp6_ = D[1*3+0];
            double tmp7_ = D[1*3+1];
            double tmp8_ = D[1*3+2];
            double tmp9_ = D[2*3+0];
            double tmp10_ = D[2*3+1];
            double tmp11_ = D[2*3+2];
            double tmp12_ = w;
            double tmp13_ = eps[0*3+0];
            double tmp14_ = eps[0*3+1];
            double tmp15_ = eps[0*3+2];
            double tmp16_ = f[0*2+0];
            double tmp17_ = f[0*2+1];
            double tmp18_ = tmp1_ * tmp3_;
            double tmp19_ = tmp2_ * tmp5_;
            double tmp20_ = tmp18_ + tmp19_;
            double tmp21_ = tmp2_ * tmp4_;
            double tmp22_ = tmp1_ * tmp5_;
            double tmp23_ = tmp21_ + tmp22_;
            double tmp24_ = tmp1_ * tmp6_;
            double tmp25_ = tmp2_ * tmp8_;
            double tmp26_ = tmp24_ + tmp25_;
            double tmp27_ = tmp2_ * tmp7_;
            double tmp28_ = tmp1_ * tmp8_;
            double tmp29_ = tmp27_ + tmp28_;
            double tmp30_ = tmp1_ * tmp9_;
            double tmp31_ = tmp2_ * tmp11_;
            double tmp32_ = tmp30_ + tmp31_;
            double tmp33_ = tmp2_ * tmp10_;
            double tmp34_ = tmp1_ * tmp11_;
            double tmp35_ = tmp33_ + tmp34_;
            double tmp36_ = tmp20_ * tmp13_;
            double tmp37_ = tmp26_ * tmp14_;
            double tmp38_ = tmp36_ + tmp37_;
            double tmp39_ = tmp32_ * tmp15_;
            double tmp40_ = tmp38_ + tmp39_;
            double tmp41_ = tmp23_ * tmp13_;
            double tmp42_ = tmp29_ * tmp14_;
            double tmp43_ = tmp41_ + tmp42_;
            double tmp44_ = tmp35_ * tmp15_;
            double tmp45_ = tmp43_ + tmp44_;
            double tmp46_ = tmp40_ * tmp12_;
            double tmp47_ = tmp45_ * tmp12_;
            double tmp48_ = tmp16_ + tmp46_;
            double tmp49_ = tmp17_ + tmp47_;
            f[0*2+0] = tmp48_;
            f[0*2+1] = tmp49_;
            } /* </generated matexpr> */
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
                /* <generated matexpr> */ {
                double tmp1_ = D[0*3+0];
                double tmp2_ = D[0*3+1];
                double tmp3_ = D[0*3+2];
                double tmp4_ = D[1*3+0];
                double tmp5_ = D[1*3+1];
                double tmp6_ = D[1*3+2];
                double tmp7_ = D[2*3+0];
                double tmp8_ = D[2*3+1];
                double tmp9_ = D[2*3+2];
                double tmp10_ = w;
                double tmp11_ = Ni_x;
                double tmp12_ = Ni_y;
                double tmp13_ = Nj_x;
                double tmp14_ = Nj_y;
                double tmp15_ = Knodal[0*8+0];
                double tmp16_ = Knodal[0*8+1];
                double tmp17_ = Knodal[1*8+0];
                double tmp18_ = Knodal[1*8+1];
                double tmp19_ = tmp11_ * tmp1_;
                double tmp20_ = tmp12_ * tmp3_;
                double tmp21_ = tmp19_ + tmp20_;
                double tmp22_ = tmp12_ * tmp2_;
                double tmp23_ = tmp11_ * tmp3_;
                double tmp24_ = tmp22_ + tmp23_;
                double tmp25_ = tmp11_ * tmp4_;
                double tmp26_ = tmp12_ * tmp6_;
                double tmp27_ = tmp25_ + tmp26_;
                double tmp28_ = tmp12_ * tmp5_;
                double tmp29_ = tmp11_ * tmp6_;
                double tmp30_ = tmp28_ + tmp29_;
                double tmp31_ = tmp11_ * tmp7_;
                double tmp32_ = tmp12_ * tmp9_;
                double tmp33_ = tmp31_ + tmp32_;
                double tmp34_ = tmp12_ * tmp8_;
                double tmp35_ = tmp11_ * tmp9_;
                double tmp36_ = tmp34_ + tmp35_;
                double tmp37_ = tmp21_ * tmp13_;
                double tmp38_ = tmp33_ * tmp14_;
                double tmp39_ = tmp37_ + tmp38_;
                double tmp40_ = tmp24_ * tmp13_;
                double tmp41_ = tmp36_ * tmp14_;
                double tmp42_ = tmp40_ + tmp41_;
                double tmp43_ = tmp27_ * tmp14_;
                double tmp44_ = tmp33_ * tmp13_;
                double tmp45_ = tmp43_ + tmp44_;
                double tmp46_ = tmp30_ * tmp14_;
                double tmp47_ = tmp36_ * tmp13_;
                double tmp48_ = tmp46_ + tmp47_;
                double tmp49_ = tmp39_ * tmp10_;
                double tmp50_ = tmp42_ * tmp10_;
                double tmp51_ = tmp45_ * tmp10_;
                double tmp52_ = tmp48_ * tmp10_;
                double tmp53_ = tmp15_ + tmp49_;
                double tmp54_ = tmp16_ + tmp50_;
                double tmp55_ = tmp17_ + tmp51_;
                double tmp56_ = tmp18_ + tmp52_;
                Knodal[0*8+0] = tmp53_;
                Knodal[0*8+1] = tmp54_;
                Knodal[1*8+0] = tmp55_;
                Knodal[1*8+1] = tmp56_;
                } /* </generated matexpr> */
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
