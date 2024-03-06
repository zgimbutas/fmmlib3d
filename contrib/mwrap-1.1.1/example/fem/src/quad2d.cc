/*
 * quad2d.h
 *   Implementation for four-node quad element shapes.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */
#include "quad2d.h"

#define ME Quad2d


ME::~ME()
{
}


void ME::set_node(int nodenum, const double* x)
{
    nodex[2*nodenum+0] = x[0];
    nodex[2*nodenum+1] = x[1];
}


void ME::set_nodes(const double* x)
{
    for (int i = 0; i < 8; ++i)
        nodex[i] = x[i];
}


void ME::eval(const double *XX, double* xx,
              double* N, double* dN,
              double& J) const
{
    double X = XX[0];
    double Y = XX[1];
    double FF[4];

    /* <generator matexpr>
    // Evaluate element functions

    input X, Y;
    input nodex(2,4);

    output N(4);
    output dN(4,2);
    output xx(2);
    output FF(2,2);

    N1x = (1-X)/2;  N1y = (1-Y)/2;
    N2x = (1+X)/2;  N2y = (1+Y)/2;    

    N  = [ N1x*N1y; N2x*N1y; N2x*N2y; N1x*N2y];
    dN = [    -N1y,     N1y,     N2y,    -N2y;
              -N1x,    -N2x,     N2x,     N1x ]'/2;

    xx = nodex*N;
    FF = nodex*dN;
    */
    remap_gradients(FF, dN, J);
}


void ME::remap_gradients(double* FF, double* dN, double& J) const
{
    J = (FF[0]*FF[3]-FF[1]*FF[2]);
    double invF[4] = {  FF[3]/J, -FF[1]/J, 
		       -FF[2]/J,  FF[0]/J };

    /* <generator matexpr>
    // Remap gradients

    inout dN(4,2);
    input invF(2,2);
    dN = dN*invF;
    */
}
