/*
 * mwrap-mgen.cc
 *   Generate MATLAB scripts from MWrap AST.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cassert>
#include "mwrap-ast.h"


/* -- Output MATLAB call code -- */
/*
 * Each call is translated to MATLAB code of the form
 *
 *   mex_id_ = 'identifier from id_string';
 *   [result1, result2, ...] = mexfunc(mex_id_, in1, ..., dim1, dim2, ...);
 *
 * where dim1, dim2, ... are any explicitly provided array dimensions.
 */


int has_output_args(Var* v)
{
    if (!v)
        return 0;
    if (v->iospec == 'o' || v->iospec == 'b')
        return 1;
    else 
        return has_output_args(v->next);
}


void print_output_args(FILE* fp, Var* v, int& first)
{
    if (!v)
        return;
    if (v->iospec == 'o' || v->iospec == 'b') {
        if (!first)
            fprintf(fp, ", ");
        fprintf(fp, "%s", v->name);
        first = 0;
    }
    print_output_args(fp, v->next, first);
}


void print_input_args(FILE* fp, Var* v)
{
    if (!v)
        return;
    if (v->tinfo == VT_const)
        fprintf(fp, ", 0", v->name);
    else if (v->iospec == 'i' || v->iospec == 'b')
        fprintf(fp, ", %s", v->name);
    print_input_args(fp, v->next);
}


void print_dimension_args(FILE* fp, Expr* e)
{
    if (!e)
        return;
    fprintf(fp, ", %s", e->value);
    print_dimension_args(fp, e->next);
}


void print_dimension_args(FILE* fp, Var* v)
{
    if (!v)
        return;
    if (v->qual)
        print_dimension_args(fp, v->qual->args);
    print_dimension_args(fp, v->next);
}


void print_matlab_call(FILE* fp, Func* f, const char* mexfunc)
{
    fprintf(fp, "mex_id_ = '%s';\n", id_string(f).c_str());
    if (f->ret || has_output_args(f->args)) {
        int first = 1;
        fprintf(fp, "[");
        print_output_args(fp, f->ret, first);
        print_output_args(fp, f->args, first);
        fprintf(fp, "] = ");
    }
    fprintf(fp, "%s(mex_id_", mexfunc);
    if (f->thisv)
        fprintf(fp, ", %s", f->thisv);
    print_input_args(fp, f->args);
    print_dimension_args(fp, f->ret);
    print_dimension_args(fp, f->args);
    fprintf(fp, ");\n");
}
