/*
 * mwrap-typecheck.cc
 *   Typecheck MWrap AST.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cassert>
#include "mwrap-ast.h"


/* -- Label input / output indices -- */
/*
 * To each argument, we assign an input label (index of the argument
 * in the prhs array) and an output label (the index in the plhs
 * array).  The protocol is that we pass input and inout arguments in
 * first, and then tack on all the dimensioning arguments at the end.
 */


void label_dim_args(Expr* e, int& icount)
{
    if (!e)
        return;
    e->input_label = icount++;
    label_dim_args(e->next, icount);
}


void label_dim_args(Var* v, int& icount)
{
    if (!v)
        return;
    if (v->qual)
        label_dim_args(v->qual->args, icount);
    label_dim_args(v->next, icount);
}


void label_args(Var* v, int& icount, int& ocount)
{
    if (!v)
        return;
    if (v->iospec == 'i' || v->iospec == 'b')
        v->input_label = icount++;
    if (v->iospec == 'o' || v->iospec == 'b')
        v->output_label = ocount++;
    label_args(v->next, icount, ocount);
}


void label_args(Func* f)
{
    int icount = 0;
    int ocount = 0;
    if (f->thisv)
        icount = 1;
    label_args(f->ret, icount, ocount);
    label_args(f->args, icount, ocount);
    label_dim_args(f->ret, icount);
    label_dim_args(f->args, icount);
}


/* -- Fortran-ize arguments -- */
/*
 * For FORTRAN calls, we make the following type conversions and
 * checks on arguments:
 *
 *  - All types of scalar arguments are converted to pointer-to-scalar.
 *  - Object arguments are forbidden
 *  - C string arguments generate a warning
 *
 * Only scalar return values are allowed from wrapped FORTRAN functions.
 */


int fortranize_args(Var* v, int line)
{
    if (!v)
        return 0;
    int errcount = 0;
    if (v->tinfo == VT_obj || v->tinfo == VT_p_obj || v->tinfo == VT_r_obj) {
        fprintf(stderr, "Error (%d): ", line);
        fprintf(stderr, "Cannot pass object %s to FORTRAN\n", v->name);
        ++errcount;
    } else if (v->tinfo == VT_rarray) {
        fprintf(stderr, "Error (%d): ", line);
        fprintf(stderr, "Cannot pass pointer ref %s to FORTRAN\n", v->name);
        ++errcount;
    } else if (v->tinfo == VT_string) {
        fprintf(stderr, "Warning (%d): ", line);
        fprintf(stderr, "Danger passing C string %s to FORTRAN\n", v->name);
    } else if (v->tinfo == VT_scalar || v->tinfo == VT_r_scalar) {
        v->tinfo = VT_p_scalar;
    } else if (v->tinfo == VT_cscalar || v->tinfo == VT_r_cscalar) {
        v->tinfo = VT_p_cscalar;
    } else if (v->tinfo == VT_zscalar || v->tinfo == VT_r_zscalar) {
        v->tinfo = VT_p_zscalar;
    }
    return errcount + fortranize_args(v->next, line);
}


int fortranize_ret(Var* v, int line)
{
    if (!v)
        return 0;
    if (v->tinfo == VT_cscalar || v->tinfo == VT_zscalar) {
        fprintf(stderr, "Warning (%d): ", line);
        fprintf(stderr, "Danger returning complex from FORTRAN\n");
    } else if (v->tinfo != VT_scalar) {
        fprintf(stderr, "Error (%d): ", line);
        fprintf(stderr, "Can only return scalars from FORTRAN\n");
        return 1;
    }
    return 0;
}


int fortranize_args(Func* f, int line)
{
    if (!f->fort)
        return 0;
    return (fortranize_args(f->args, line) +
            fortranize_ret(f->ret, line));
}


/* -- Type info assignment and checking -- */
/*
 * Strictly speaking, we're doing semantic checking here -- there are a
 * few things we check for that don't qualify as type errors.  The
 * assumption is that any input to the C code generator has passed the
 * type checker.
 */


int assign_scalar_tinfo(Var* v, int line,
                        int tags, int tagp, int tagr, int taga,
                        int tagar)
{
    // Numeric types
    if (!v->qual) {
        v->tinfo = tags;
    } else if (v->qual->qual == '*') {
        v->tinfo = tagp;
    } else if (v->qual->qual == '&') {
        v->tinfo = tagr;
    } else if (v->qual->qual == 'a') {
        v->tinfo = taga;
        if (v->qual->args && 
            v->qual->args->next && 
            v->qual->args->next->next) {
            fprintf(stderr, "Error (%d): ", line);
            fprintf(stderr, "Array %s should be 1D or 2D\n", v->name);
            return 1;
        }
    } else if (v->qual->qual == 'r') {
        v->tinfo = tagar;
        if (tagar == VT_unk) {
            fprintf(stderr, "Error (%d): ", line);
            fprintf(stderr, "Array ref %s must be to a real array\n", 
                    v->name);
            return 1;
        }
        if (v->qual->args && 
            v->qual->args->next && 
            v->qual->args->next->next) {
            fprintf(stderr, "Error (%d): ", line);
            fprintf(stderr, "Array %s should be 1D or 2D\n", v->name);
            return 1;
        }
    } else {
        assert(0);
    }

    return 0;
}


int assign_tinfo(Var* v, int line)
{
    if (is_scalar_type(v->basetype)) {

        return assign_scalar_tinfo(v, line, VT_scalar, VT_p_scalar, 
                                   VT_r_scalar, VT_array, VT_rarray);

    } else if (is_cscalar_type(v->basetype)) {

        return assign_scalar_tinfo(v, line, VT_cscalar, VT_p_cscalar, 
                                   VT_r_cscalar, VT_carray, VT_unk);

    } else if (is_zscalar_type(v->basetype)) {

        return assign_scalar_tinfo(v, line, VT_zscalar, VT_p_zscalar, 
                                   VT_r_zscalar, VT_zarray, VT_unk);

    } else if (strcmp(v->basetype, "const") == 0) {

        if (v->qual) {
            fprintf(stderr, "Error (%d): ", line);
            fprintf(stderr, "Constant %s cannot have modifiers\n", v->name);
            return 1;
        }
        v->tinfo = VT_const;
        if (v->name[0] == '\'') {
            char* p_out = v->name;
            char* p_in  = v->name;
            for (; *p_in; ++p_in) {
                if (*p_in != '\'')
                    *p_out++ = *p_in;
            }
            *p_out++ = 0;
        }

    } else if (strcmp(v->basetype, "cstring") == 0) {

        // String type
        if (v->qual && v->qual->qual != 'a') {
            fprintf(stderr, "Error (%d): ", line);
            fprintf(stderr, "String type %s cannot have modifiers\n",
                    v->name);
            return 1;
        }
        if (v->qual && v->qual->args && v->qual->args->next) {
            fprintf(stderr, "Error (%d): ", line);
            fprintf(stderr, "Strings are one dimensional\n", v->name);
            return 1;
        }
        v->tinfo = VT_string;

    } else if (strcmp(v->basetype, "mxArray") == 0) {

        // MATLAB intrinsic type
        if (v->qual) {
            fprintf(stderr, "Error (%d): ", line);
            fprintf(stderr, "mxArray %s cannot have modifiers\n",
                    v->name);
            return 1;
        }
        v->tinfo = VT_mx;

    } else {

        // Object type
        if (!v->qual)
            v->tinfo = VT_obj;
        else if (v->qual->qual == '*')
            v->tinfo = VT_p_obj;
        else if (v->qual->qual == '&')
            v->tinfo = VT_r_obj;
        else if (v->qual->qual == 'a' || v->qual->qual == 'r') {
            fprintf(stderr, "Error (%d): ", line);
            fprintf(stderr, "%s cannot be an array of object %s\n",
                    v->name, v->basetype);
            return 1;
        } else
            assert(0);

    }
    return 0;
}


int typecheck_return(Var* v, int line)
{
    if (!v)
        return 0;

    int err = assign_tinfo(v, line);
    if ((v->tinfo == VT_array ||
         v->tinfo == VT_carray ||
         v->tinfo == VT_zarray) && 
        !(v->qual && v->qual->args)) {
        fprintf(stderr, "Error (%d): ", line);
        fprintf(stderr, "Return array %s must have dims\n", v->name);
        ++err;
    } else if (v->tinfo == VT_const) {
        fprintf(stderr, "Error (%d): ", line);
        fprintf(stderr, "Cannot return constant\n");
        ++err;
    } else if (v->tinfo == VT_rarray) {
        fprintf(stderr, "Error (%d): ", line);
        fprintf(stderr, "Ref to array %s looks just like array on return\n", 
                v->name);
        ++err;
    } else if (v->tinfo == VT_string && v->qual) {
        fprintf(stderr, "Error (%d): ", line);
        fprintf(stderr, "Return string %s cannot have dims\n", v->name);
        ++err;
    }
    return err;
}


int typecheck_args(Var* v, int line)
{
    if (!v)
        return 0;

    int err = assign_tinfo(v, line);
    if (v->iospec == 'i')
        return err + typecheck_args(v->next, line);

    if (isdigit(v->name[0])) {
        fprintf(stderr, "Error (%d): ", line);
        fprintf(stderr, "Number %s cannot be output\n", v->name);
        ++err;
    }

    if ((v->tinfo == VT_obj || v->tinfo == VT_p_obj || v->tinfo == VT_r_obj) &&
        !is_mxarray_type(v->basetype)) {
        fprintf(stderr, "Error (%d): ", line);
        fprintf(stderr, "Object %s cannot be output\n", v->name);
        ++err;
    } else if ((v->tinfo == VT_array || 
                v->tinfo == VT_carray ||
                v->tinfo == VT_zarray ||
                v->tinfo == VT_rarray) && 
               v->iospec == 'o' && 
               !(v->qual && v->qual->args)) {
        fprintf(stderr, "Error (%d): ", line);
        fprintf(stderr, "Output array %s must have dims\n", v->name);
        ++err;
    } else if (v->tinfo == VT_rarray && v->iospec != 'o') {
        fprintf(stderr, "Error (%d): ", line);
        fprintf(stderr, "Array ref %s *must* be output\n", v->name);
        ++err;
    } else if (v->tinfo == VT_scalar) {
        fprintf(stderr, "Error (%d): ", line);
        fprintf(stderr, "Scalar %s cannot be output\n", v->name);
        ++err;
    } else if (v->tinfo == VT_const) {
        fprintf(stderr, "Error (%d): ", line);
        fprintf(stderr, "Constant %s cannot be output\n", v->name);
        ++err;
    } else if (v->tinfo == VT_string && 
               !(v->qual && v->qual->args)) {
        fprintf(stderr, "Error (%d): ", line);
        fprintf(stderr, "String %s cannot be output without size\n", v->name);
        ++err;
    } else if (v->tinfo == VT_mx && v->iospec == 'b') {
        fprintf(stderr, "Error (%d): ", line);
        fprintf(stderr, "mxArray %s cannot be used for inout\n", v->name);
        ++err;
    }
    return err + typecheck_args(v->next, line);
}


int typecheck(Func* f, int line)
{
    label_args(f);
    return (typecheck_return(f->ret, line) + 
            typecheck_args(f->args, line) +
            fortranize_args(f, line));
}

