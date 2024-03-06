/*
 * quad2d.h
 *   Interface for four-node quad element shapes.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */
#ifndef QUAD2D_H
#define QUAD2D_H

#include "feshapes.h"

class Quad2d : public FEShapes {
public:
    Quad2d() {}
    Quad2d(const double* x) { set_nodes(x); }
    virtual ~Quad2d();

    void set_node(int nodenum, const double* x);
    void set_nodes(const double* x);

    void eval(const double *XX, double* xx,
              double* N, double* dN,
              double& J) const;

    int nshape() const { return 4; }

private:
    double nodex[2*4];
    void remap_gradients(double* FF, double* dN, double& J) const;
};

#endif /* QUAD2D_H */
