/*
 * gauss2by2.cc
 *   Parent domain node locations and weights for quadrature on a
 *   2-by-2 Gauss grid over [-1,1]^2.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include "gauss2by2.h"

const double Gauss4::X[12] = {
    -0.577350269189626, -0.577350269189626, 1,
     0.577350269189626, -0.577350269189626, 1,
     0.577350269189626,  0.577350269189626, 1,
    -0.577350269189626,  0.577350269189626, 1,
};
