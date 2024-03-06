/*
 * feshapes.h
 *   Interface for element shape functions.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#ifndef FESHAPES_H
#define FESHAPES_H

class FEShapes {
public:
    virtual ~FEShapes() {}

    virtual void set_node(int nodenum, const double* x) = 0;
    virtual void set_nodes(const double* x) = 0;

    virtual void eval(const double *XX, double* xx,
                      double* N, double* dN,
                      double& J) const = 0;

    virtual int nshape() const = 0;
};

#endif /* FESHAPES_H */
