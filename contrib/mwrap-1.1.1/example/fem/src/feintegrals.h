/*
 * feintegrals.h
 *   Interface for element quadrature rules.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#ifndef FEINTEGRALS_H
#define FEINTEGRALS_H

#include "feshapes.h"

#include <vector>
using std::vector;

class VolumeQuadrature {
public:
    VolumeQuadrature(FEShapes& quad, int ndim) : 
        quad(quad), nshape1(quad.nshape()),
        xx1(2), N1(quad.nshape()), dN1(ndim*quad.nshape()) {}

    virtual void start() = 0;
    virtual bool done() = 0;

    virtual void operator++() = 0;
    virtual double wt() = 0;

    int nshape()            { return nshape1;   }
    double* xx()            { return &(xx1[0]); }
    double* N()             { return &(N1[0]);  }
    double* dN()            { return &(dN1[0]); }
    double xx(int i)        { return xx1[i]; }
    double N(int i)         { return N1[i];   }
    double dN(int i, int j) { return dN1[i+j*nshape1]; }

protected:
    FEShapes& quad;
    int nshape1;

    vector<double> xx1;
    vector<double> N1;
    vector<double> dN1;
    double J;
};

#endif /* FEINTEGRALS_H */
