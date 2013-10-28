/*
 * gauss2by2.h
 *   2-by-2 Gauss grid quadrature over [-1,1]^2.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#ifndef GAUSS2BY2_H
#define GAUSS2BY2_H

#include "feshapes.h"
#include "feintegrals.h"

class Gauss4 : public VolumeQuadrature {
public:
    Gauss4(FEShapes& quad) : VolumeQuadrature(quad, 2) {
        start();
    }

    void start()      { i = 0; eval();        }
    bool done()       { return (i > 3);       }
    void operator++() { if (++i <= 3) eval(); }
    double wt()       { return X[3*i+2]*J;    }

private:
    static const double X[12]; 
    int i;

    void eval() { quad.eval(X+3*i, xx(), N(), dN(), J); }
};

#endif /* GAUSS2BY2_H */
