/*
 * etype.h
 *   Element type interface definition.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#ifndef ETYPE_H
#define ETYPE_H

#include "assembler.h"

class Mesh;

class EType {
public:
    virtual ~EType();
    virtual void assign_ids(Mesh* mesh, int eltid) = 0;
    virtual void assemble_f(Mesh* mesh, int eltid) = 0;
    virtual void assemble_K(Mesh* mesh, int eltid, 
                            MatrixAssembler* K_assembler) = 0;
};

#endif /* ETYPE_H */
