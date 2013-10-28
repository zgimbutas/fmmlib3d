/*
 * mwrap-cgen.cc
 *   Generate C MEX file from MWrap AST.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cassert>
#include "mwrap-ast.h"
#include "mwrap-support.h"


/* -- General utility functions -- */

/*
 * Check whether a tinfo is an array type
 */
bool is_array(int tinfo)
{
    return (tinfo == VT_array || tinfo == VT_carray || tinfo == VT_zarray);
}


/*
 * Check whether a tinfo is an array type
 */
bool is_obj(int tinfo)
{
    return (tinfo == VT_obj || tinfo == VT_p_obj || tinfo == VT_r_obj);
}


/*
 * A function has a nullable return value if it returns a non-object
 * type for which a NULL pointer value makes sense.
 */
bool nullable_return(Func* f)
{
    return (f->ret &&
            (f->ret->tinfo == VT_string ||
             is_array(f->ret->tinfo) ||
             f->ret->tinfo == VT_p_scalar ||
             f->ret->tinfo == VT_p_cscalar ||
             f->ret->tinfo == VT_p_zscalar));
}


/*
 * Get the name of a variable into the provided buffer.
 */
const char* vname(Var* v, char* buf)
{
    if (v->iospec == 'o')
        sprintf(buf, "out%d_", v->output_label);
    else
        sprintf(buf, "in%d_", v->input_label);
    return buf;
}


/*
 * Check whether there are any FORTRAN routines in the list.
 */
bool has_fortran(Func* f)
{
    for (; f; f=f->next)
        if (f->fort)
            return true;
    return false;
}


/*
 * Check whether the specified variable is complex
 */
bool complex_tinfo(Var* v)
{
    return (v->tinfo == VT_carray || v->tinfo == VT_zarray ||
            v->tinfo == VT_cscalar || v->tinfo == VT_zscalar ||
            v->tinfo == VT_r_cscalar || v->tinfo == VT_r_zscalar ||
            v->tinfo == VT_p_cscalar || v->tinfo == VT_p_zscalar);
}


/*
 * Count routines in the list.
 */
int max_routine_id(Func* f)
{
    int maxid = 0;
    for (; f; f=f->next)
        if (f->id > maxid)
            maxid = f->id;
    return maxid;
}


/* -- Generate standard complex type definitions -- */

void mex_cpp_complex(FILE* fp)
{
    fprintf(fp, 
            "#include <complex>\n"
            "\n"
            "typedef std::complex<double> dcomplex;\n"
            "#define real_dcomplex(z) std::real(z)\n"
            "#define imag_dcomplex(z) std::imag(z)\n"
            "#define setz_dcomplex(z,r,i)  *z = dcomplex(r,i)\n"
            "\n"
            "typedef std::complex<float> fcomplex;\n"
            "#define real_fcomplex(z) std::real(z)\n"
            "#define imag_fcomplex(z) std::imag(z)\n"
            "#define setz_fcomplex(z,r,i)  *z = fcomplex(r,i)\n\n");
}


void mex_c99_complex(FILE* fp)
{
    fprintf(fp, 
            "#include <complex.h>\n"
            "\n"
            "typedef _Complex double dcomplex;\n"
            "#define real_dcomplex(z) creal(z)\n"
            "#define imag_dcomplex(z) cimag(z)\n"
            "#define setz_dcomplex(z,r,i)  *z = r + i*_Complex_I\n"
            "\n"
            "typedef _Complex float fcomplex;\n"
            "#define real_fcomplex(z) crealf(z)\n"
            "#define imag_fcomplex(z) cimagf(z)\n"
            "#define setz_fcomplex(z,r,i)  *z = r + i*_Complex_I\n\n");
}


/* -- Generate scalar conversion info -- */
/*
 * For each scalar type T, we define methods to convert a double array
 * to/from a T array, and to handle return of a pointer to a single T
 * (which might be NULL).
 */


void mex_define_copiers(FILE* fp, const char* name)
{
    fprintf(fp, "mxWrapGetArrayDef(mxWrapGetArray_%s, %s)\n", name, name);
    fprintf(fp, "mxWrapCopyDef    (mxWrapCopy_%s,     %s)\n", name, name);
    fprintf(fp, "mxWrapReturnDef  (mxWrapReturn_%s,   %s)\n", name, name);
}


void mex_define_zcopiers(FILE* fp, const char* name, const char* ztype)
{
    fprintf(fp, 
            "mxWrapGetScalarZDef(mxWrapGetScalar_%s, %s,\n"
            "                    %s, setz_%s)\n", 
            name, name, ztype, name);
    fprintf(fp, 
            "mxWrapGetArrayZDef (mxWrapGetArray_%s, %s,\n"
            "                    %s, setz_%s)\n", 
            name, name, ztype, name);
    fprintf(fp, 
            "mxWrapCopyZDef     (mxWrapCopy_%s, %s,\n"
            "                    real_%s, imag_%s)\n", 
            name, name, name, name);
    fprintf(fp, 
            "mxWrapReturnZDef   (mxWrapReturn_%s, %s,\n"
            "                    real_%s, imag_%s)\n", 
            name, name, name, name);
}


void mex_define_copiers(FILE* fp)
{
    fprintf(fp, "/* Array copier definitions */\n");
    for( set<string>::iterator iter = scalar_decls.begin();
         iter != scalar_decls.end(); 
         ++iter)
        mex_define_copiers(fp, iter->c_str());
    for( set<string>::iterator iter = cscalar_decls.begin();
         iter != cscalar_decls.end(); 
         ++iter)
        mex_define_zcopiers(fp, iter->c_str(), "float");
    for( set<string>::iterator iter = zscalar_decls.begin();
         iter != zscalar_decls.end(); 
         ++iter)
        mex_define_zcopiers(fp, iter->c_str(), "double");
    fprintf(fp, "\n");
}


/* -- Handle FORTRAN name mangling -- */
/*
 * For each FORTRAN function, we define a macro version of the name
 * starting with MWF77_, which will handle the various name-mangling
 * issues that come up.  We support three name mangling conventions for
 * now: all upper (no underscore mangling); all lower, single underscore;
 * and all lower, single underscore, double underscore for names that
 * contain an underscore (f2c).
 */


void mex_define_caps_fname(FILE* fp, Func* f)
{
    fprintf(fp, "#define MWF77_%s ", f->funcv);
    for (char* s = f->funcv; *s; ++s)
        fputc(toupper(*s), fp);
    fprintf(fp, "\n");
}


void mex_define_caps_fnames(FILE* fp, Func* f)
{
    set<string> fnames;
    for (; f; f = f->next) {
        if (f->fort && fnames.find(f->funcv) == fnames.end()) {
            mex_define_caps_fname(fp, f);
            fnames.insert(f->funcv);
        }
    }
}


void mex_define_underscore1_fname(FILE* fp, Func* f)
{
    fprintf(fp, "#define MWF77_%s ", f->funcv);
    for (char* s = f->funcv; *s; ++s)
        fputc(tolower(*s), fp);
    fprintf(fp, "_\n");
}


void mex_define_underscore1_fnames(FILE* fp, Func* f)
{
    set<string> fnames;
    for (; f; f = f->next) {
        if (f->fort && fnames.find(f->funcv) == fnames.end()) {
            mex_define_underscore1_fname(fp, f);
            fnames.insert(f->funcv);
        }
    }
}


void mex_define_f2c_fname(FILE* fp, Func* f)
{
    fprintf(fp, "#define MWF77_%s ", f->funcv);
    bool has_underscore = false;
    for (char* s = f->funcv; *s; ++s) {
        has_underscore = (has_underscore || (*s == '_'));
        fputc(tolower(*s), fp);
    }
    if (has_underscore)
        fprintf(fp, "__\n");
    else
        fprintf(fp, "_\n");
}


void mex_define_f2c_fnames(FILE* fp, Func* f)
{
    set<string> fnames;
    for (; f; f = f->next) {
        if (f->fort && fnames.find(f->funcv) == fnames.end()) {
            mex_define_f2c_fname(fp, f);
            fnames.insert(f->funcv);
        }
    }
}


void mex_define_fnames(FILE* fp, Func* f)
{
    fprintf(fp, "#if defined(MWF77_CAPS)\n");
    mex_define_caps_fnames(fp, f);
    fprintf(fp, "#elif defined(MWF77_UNDERSCORE1)\n");
    mex_define_underscore1_fnames(fp, f);
    fprintf(fp, "#else /* f2c convention */\n");
    mex_define_f2c_fnames(fp, f);
    fprintf(fp, "#endif\n\n");
}


/* -- Handle FORTRAN declarations -- */
/*
 * We assume unless told otherwise that all FORTRAN functions need
 * prototypes, which we will generate.
 */


void mex_fortran_arg(FILE* fp, Var* v)
{
    if (!v)
        return;
    if (v->tinfo == VT_mx)
        fprintf(fp, "const mxArray*", v->name);
    else
        fprintf(fp, "%s*", v->basetype, v->name);
    if (v->next) {
        fprintf(fp, ", ");
        mex_fortran_arg(fp, v->next);
    }
}


void mex_fortran_decl(FILE* fp, Func* f)
{
    if (f->ret)
        fprintf(fp, "%s ", f->ret->basetype);
    else
        fprintf(fp, "MWF77_RETURN ");
    fprintf(fp, "MWF77_%s(", f->funcv);
    mex_fortran_arg(fp, f->args);
    fprintf(fp, ");\n");
}


void mex_fortran_decls(FILE* fp, Func* f)
{
    fprintf(fp, 
            "#ifdef __cplusplus\n"
            "extern \"C\" { /* Prevent C++ name mangling */\n"
            "#endif\n\n"
            "#ifndef MWF77_RETURN\n"
            "#define MWF77_RETURN int\n"
            "#endif\n\n");

    set<string> fnames;
    for (; f; f = f->next) {
        if (f->fort && fnames.find(f->funcv) == fnames.end()) {
            mex_fortran_decl(fp, f);
            fnames.insert(f->funcv);
        }
    }

    fprintf(fp,
            "\n#ifdef __cplusplus\n"
            "} /* end extern C */\n"
            "#endif\n\n");
}


/* -- Generate class conversion info -- */
/*
 * For every parent class, we define a getter method that can take a
 * string of the form "T:value" and interpret the value as a pointer
 * to T for every child class T.  We want to get an *exact* match with
 * the child pointer type and then cast that pointer to the parent type
 * in the C++ code, because otherwise we don't get the right behavior
 * with multiple inheritance (i.e. when a cast actually changes the pointer
 * value).
 */


void mex_casting_getter_type(FILE* fp, const char* name)
{
    fprintf(fp, 
            "    %s* p_%s = NULL;\n"
            "    sscanf(pbuf, \"%s:%%p\", &p_%s);\n"
            "    if (p_%s)\n"
            "        return p_%s;\n\n", 
            name, name, name, name, name, name);
}


void mex_casting_getter(FILE* fp, const char* cname, 
                        InheritsDecl* inherits)
{
    fprintf(fp, "\n%s* mxWrapGetP_%s(const mxArray* a, const char** e)\n", 
            cname, cname);
    fprintf(fp, 
            "{\n"
            "    char pbuf[128];\n"
            "    if (mxGetClassID(a) == mxDOUBLE_CLASS &&\n"
            "        mxGetM(a)*mxGetN(a) == 1 && *mxGetPr(a) == 0)\n"
            "        return NULL;\n"
            "    if (!mxIsChar(a)) {\n"
            "#ifdef R2008OO\n"
            "        mxArray* ap = mxGetProperty(a, 0, \"mwptr\");\n"
            "        if (ap)\n"
            "            return mxWrapGetP_%s(ap, e);\n"
            "#endif\n"
            "        *e = \"Invalid pointer\";\n"
            "        return NULL;\n"
            "    }\n"
            "    mxGetString(a, pbuf, sizeof(pbuf));\n\n", cname);

    mex_casting_getter_type(fp, cname);
    for (InheritsDecl* i = inherits; i; i = i->next)
        mex_casting_getter_type(fp, i->name);

    fprintf(fp, 
            "    *e = \"Invalid pointer to %s\";\n"
            "    return NULL;\n"
            "}\n\n", cname);
}


void mex_casting_getters(FILE* fp)
{
    map<string,InheritsDecl*>::iterator i = class_decls.begin();
    for (; i != class_decls.end(); ++i)
        mex_casting_getter(fp, i->first.c_str(), i->second);
}


void mex_cast_get_p(FILE* fp, const char* basetype, int input_label)
{
    fprintf(fp, "    in%d_ = ", input_label);
    if (is_mxarray_type(basetype))
        fprintf(fp, "mxWrapGet_%s(prhs[%d], &mw_err_txt_);\n",
                basetype, input_label);
    else if (class_decls.find(basetype) == class_decls.end())
        fprintf(fp, "(%s*) mxWrapGetP(prhs[%d], \"%s:%%p\", &mw_err_txt_);\n",
                basetype, input_label, basetype);
    else
        fprintf(fp, "mxWrapGetP_%s(prhs[%d], &mw_err_txt_);\n", 
                basetype, input_label);
    fprintf(fp, 
            "    if (mw_err_txt_)\n"
            "        goto mw_err_label;\n");
}


/* -- Mex stub variable declarations -- */
/*
 * In each stub function, we declare local variables corresponding to
 * arguments (input and output), return values, and dimension parameters. 
 * We use the names inX_, outX_, and dimX_ to represent input/inout, output
 * and dimension parameters, where X is the index in the prhs array or
 * plhs array.
 */


void mex_declare_type(char* typebuf, Var* v)
{
    if (is_obj(v->tinfo) || is_array(v->tinfo))
        sprintf(typebuf, "%s*", v->basetype);
    else if (v->tinfo == VT_rarray)
        sprintf(typebuf, "const %s*", v->basetype);
    else if (v->tinfo == VT_scalar ||
             v->tinfo == VT_cscalar ||
             v->tinfo == VT_zscalar ||
             v->tinfo == VT_r_scalar ||
             v->tinfo == VT_r_cscalar ||
             v->tinfo == VT_r_zscalar ||
             v->tinfo == VT_p_scalar ||
             v->tinfo == VT_p_cscalar ||
             v->tinfo == VT_p_zscalar)
        sprintf(typebuf, "%s", v->basetype);
    else if (v->tinfo == VT_string)
        sprintf(typebuf, "char*");
    else if (v->tinfo == VT_mx && v->iospec == 'i')
        sprintf(typebuf, "const mxArray*");
    else if (v->tinfo == VT_mx && v->iospec == 'o')
        sprintf(typebuf, "mxArray*");
    else {
        fprintf(stderr, "v->tinfo == %d; v->name = %s\n", v->tinfo, v->name);
        assert(0);
    }
}


void mex_declare_in_args(FILE* fp, Var* v)
{
    if (!v)
        return;
    if ((v->iospec == 'i' || v->iospec == 'b') && v->tinfo != VT_const) {
        char typebuf[128];
        mex_declare_type(typebuf, v);
        if (is_array(v->tinfo) || is_obj(v->tinfo) || v->tinfo == VT_string) {
            fprintf(fp, "    %-10s  in%d_ =0; /* %-10s */\n", 
                    typebuf, v->input_label, v->name);
        } else {
            fprintf(fp, "    %-10s  in%d_;    /* %-10s */\n", 
                    typebuf, v->input_label, v->name);
        }
    }
    mex_declare_in_args(fp, v->next);
}


void mex_declare_out_args(FILE* fp, Var* v)
{
    if (!v)
        return;
    if (v->iospec == 'o' && v->tinfo != VT_mx) {
        char typebuf[128];
        mex_declare_type(typebuf, v);
        if (is_array(v->tinfo) || is_obj(v->tinfo) || v->tinfo == VT_string) {
            fprintf(fp, "    %-10s  out%d_=0; /* %-10s */\n", 
                    typebuf, v->output_label, v->name);
        } else {
            fprintf(fp, "    %-10s  out%d_;   /* %-10s */\n", 
                    typebuf, v->output_label, v->name);
        }
    }
    mex_declare_out_args(fp, v->next);
}


void mex_declare_dim_args(FILE* fp, Expr* e)
{
    if (!e)
        return;
    fprintf(fp, "    %-10s  dim%d_;   /* %-10s */\n", 
            "mwSize", e->input_label, e->value);
    mex_declare_dim_args(fp, e->next);
}


void mex_declare_dim_args(FILE* fp, Var* v)
{
    if (!v)
        return;
    if (v->qual)
        mex_declare_dim_args(fp, v->qual->args);
    mex_declare_dim_args(fp, v->next);
}


void mex_declare_return(FILE* fp, Var* v)
{
    mex_declare_out_args(fp, v);
}


void mex_declare_args(FILE* fp, Func* f)
{
    if (f->thisv) {
        char typebuf[128];
        strcpy(typebuf, f->classv);
        strcat(typebuf, "*");
        fprintf(fp, "    %-10s  in%d_ =0; /* %-10s */\n", 
                typebuf, 0, f->thisv);
    }
    mex_declare_in_args(fp, f->args);
    if (!nullable_return(f))
        mex_declare_return(fp, f->ret);
    mex_declare_out_args(fp, f->args);
    mex_declare_dim_args(fp, f->ret);
    mex_declare_dim_args(fp, f->args);
    if (f->ret || f->args || f->thisv)
        fprintf(fp, "\n");
}


/* -- Mex stub dimension retrieval -- */
/*
 * After declaring local variables, we copy in the dimension arguments.
 */


int mex_unpack_dims(FILE* fp, Expr* e)
{
    if (!e)
        return 0;
    fprintf(fp, 
            "    dim%d_ = (mwSize) mxWrapGetScalar(prhs[%d], &mw_err_txt_);\n",
            e->input_label, e->input_label);
    return mex_unpack_dims(fp, e->next)+1;
}


int mex_unpack_dims(FILE* fp, Var* v)
{
    if (!v)
        return 0;
    if (v->qual)
        return 
            mex_unpack_dims(fp, v->qual->args) + 
            mex_unpack_dims(fp, v->next);
    return mex_unpack_dims(fp, v->next);
}


void mex_unpack_dims(FILE* fp, Func* f)
{
    if (mex_unpack_dims(fp, f->ret) || mex_unpack_dims(fp, f->args))
        fprintf(fp, "\n");
}


/* -- Mex stub dimension checks -- */
/*
 * For input and inout arrays where the dimensions are given explicitly
 * in the call, we check that the dimensions of the MATLAB array passed
 * in agree with the dimensions in the dimension arguments.
 */


void mex_check_dims(FILE* fp, Var* v)
{
    if (!v)
        return;
    if (v->iospec != 'o' && is_array(v->tinfo) && v->qual && v->qual->args) {
        Expr* a = v->qual->args;
        if (a->next) {
            fprintf(fp,
                    "    if (mxGetM(prhs[%d]) != dim%d_ ||\n"
                    "        mxGetN(prhs[%d]) != dim%d_) {\n"
                    "        mw_err_txt_ = \"Bad argument size: %s\";\n"
                    "        goto mw_err_label;\n"
                    "    }\n\n",
                    v->input_label, a->input_label,
                    v->input_label, a->next->input_label,
                    v->name);
        } else {
            fprintf(fp,
                    "    if (mxGetM(prhs[%d])*mxGetN(prhs[%d]) != dim%d_) {\n"
                    "        mw_err_txt_ = \"Bad argument size: %s\";"
                    "        goto mw_err_label;\n"
                    "    }\n\n",
                    v->input_label, v->input_label, a->input_label,
                    v->name);
        }
    }
    mex_check_dims(fp, v->next);
}


void mex_check_dims(FILE* fp, Func* f)
{
    if (!fp)
        return;
    mex_check_dims(fp, f->args);
}


/* -- Output size data -- */
/*
 * Generate code to compute the size of an array from dim arguments.
 */


void mex_alloc_size1(FILE* fp, Expr* e)
{
    if (!e)
        return;
    fprintf(fp, "dim%d_", e->input_label);
    if (e->next) {
        fprintf(fp, "*");
        mex_alloc_size1(fp, e->next);
    }
}


void mex_alloc_size(FILE* fp, Expr* e)
{
    if (!e)
        fprintf(fp, "1");
    else
        mex_alloc_size1(fp, e);
}


/* -- Unpack input data -- */
/*
 * Generate code to copy data from the prhs array into the inX_
 * local variables in a stub function.  Note that for input strings,
 * we take any explicit dimension information as what the user wants --
 * if the string is too long to fit in that explicit dimension, our
 * generated code will return with an error message.
 */


void mex_unpack_input_array(FILE* fp, Var* v)
{
    fprintf(fp, "    if (mxGetM(prhs[%d])*mxGetN(prhs[%d]) != 0) {\n",
            v->input_label, v->input_label);
    if (strcmp(v->basetype, "double") == 0 && v->iospec == 'i')
        fprintf(fp, "        in%d_ = mxGetPr(prhs[%d]);\n",
                v->input_label, v->input_label);
    else
        fprintf(fp, 
                "        in%d_ = mxWrapGetArray_%s(prhs[%d], &mw_err_txt_);\n"
                "        if (mw_err_txt_)\n"
                "            goto mw_err_label;\n",
                v->input_label, v->basetype, v->input_label);
    fprintf(fp, 
            "    } else\n"
            "        in%d_ = NULL;\n", 
            v->input_label, v->basetype);
}


void mex_unpack_input_string(FILE* fp, Var* v)
{
    if (!(v->qual && v->qual->args))
        fprintf(fp, 
                "    in%d_ = mxWrapGetString(prhs[%d], &mw_err_txt_);\n"
                "    if (mw_err_txt_)\n"
                "        goto mw_err_label;\n",
                v->input_label, v->input_label);
    else {
        fprintf(fp, "    in%d_ = (char*) mxMalloc(",
                v->input_label);
        mex_alloc_size(fp, v->qual->args);
        fprintf(fp, "*sizeof(char));\n");
        fprintf(fp, "    if (mxGetString(prhs[%d], in%d_, ",
                v->input_label, v->input_label);
        mex_alloc_size(fp, v->qual->args);
        fprintf(fp, ") != 0) {\n"
                "        mw_err_txt_ = \"Invalid string argument\";\n"
                "        goto mw_err_label;\n"
                "    }\n");
    }
}


void mex_unpack_inputs(FILE* fp, Var* v)
{
    if (!v)
        return;
    if (v->iospec == 'o') {
        mex_unpack_inputs(fp, v->next);
        return;
    }

    if (is_obj(v->tinfo))
        mex_cast_get_p(fp, v->basetype, v->input_label);
    else if (is_array(v->tinfo))
        mex_unpack_input_array(fp, v);
    else if (v->tinfo == VT_scalar ||
             v->tinfo == VT_r_scalar || 
             v->tinfo == VT_p_scalar)
        fprintf(fp, 
                "    in%d_ = (%s) mxWrapGetScalar(prhs[%d], &mw_err_txt_);\n"
                "    if (mw_err_txt_)\n"
                "        goto mw_err_label;\n",
                v->input_label, v->basetype, v->input_label);
    else if (v->tinfo == VT_cscalar   || v->tinfo == VT_zscalar   ||
             v->tinfo == VT_r_cscalar || v->tinfo == VT_r_zscalar ||
             v->tinfo == VT_p_cscalar || v->tinfo == VT_p_zscalar)
        fprintf(fp, "    mxWrapGetScalar_%s(&in%d_, prhs[%d]);\n",
                v->basetype, v->input_label, v->input_label);
    else if (v->tinfo == VT_string)
        mex_unpack_input_string(fp, v);
    else if (v->tinfo == VT_mx)
        fprintf(fp, "    in%d_ = prhs[%d];\n",
                v->input_label, v->input_label);

    mex_unpack_inputs(fp, v->next);
}


void mex_unpack_inputs(FILE* fp, Func* f)
{
    if (f->thisv)
        mex_cast_get_p(fp, f->classv, 0);
    mex_unpack_inputs(fp, f->args);
}


/* -- Check input data -- */
/*
 * If an input variable is declared to be an object or a reference to
 * an object, we check to make sure the pointer passed in is not NULL.
 */


void mex_check_inputs(FILE* fp, Var* v)
{
    if (!v)
        return;
    if (v->iospec != 'o' && (v->tinfo == VT_obj || v->tinfo == VT_r_obj))
        fprintf(fp, 
                "    if (!in%d_) {\n"
                "        mw_err_txt_ = \"Argument %s cannot be null\";\n"
                "        goto mw_err_label;\n"
                "    }\n",
                v->input_label, v->name);
    mex_check_inputs(fp, v->next);
}


void mex_check_inputs(FILE* fp, Func* f)
{
    mex_check_inputs(fp, f->args);
}


/* -- Allocate output data -- */
/*
 * When we allocate data for output variables, we only allocate data
 * for *pure* output variables.  For input/output variables, we'll use
 * the scratch copy of the input as the scratch copy of the output, too.
 */


void mex_alloc_output(FILE* fp, Var* v, bool return_flag)
{
    if (!v)
        return;
    if (v->iospec == 'o') {
        if (!return_flag && is_obj(v->tinfo) && is_mxarray_type(v->basetype)) {
            fprintf(fp, "    out%d_ = mxWrapAlloc_%s();\n",
                    v->output_label, v->basetype);
        } else if (is_array(v->tinfo)) {
            fprintf(fp, "    out%d_ = (%s*) mxMalloc(",
                    v->output_label, v->basetype);
            mex_alloc_size(fp, v->qual->args);
            fprintf(fp, "*sizeof(%s));\n", v->basetype);
        } else if (v->tinfo == VT_rarray) {
            fprintf(fp, "    out%d_ = (%s*) NULL;\n",
                    v->output_label, v->basetype);
        } else if (v->tinfo == VT_string) {
            fprintf(fp, "    out%d_ = (char*) mxMalloc(", v->output_label);
            mex_alloc_size(fp, v->qual->args);
            fprintf(fp, "*sizeof(char));\n");
        }
    }
    mex_alloc_output(fp, v->next, return_flag);
}


void mex_alloc_output(FILE* fp, Func* f)
{
    if (!nullable_return(f))
        mex_alloc_output(fp, f->ret, true);
    mex_alloc_output(fp, f->args, false);
}


/* -- Record that this routine was called -- */
/*
 * For coverage testing, we keep a counter of the number of times each
 * stub is called.
 */

void mex_record_call(FILE* fp, Func* f)
{
    fprintf(fp, "    if (mexprofrecord_)\n"
                "        mexprofrecord_[%d]++;\n", f->id);
}


/* -- Make the call -- */
/*
 * The only non-obvious aspect of the call line is the way that return
 * values are treated.  The biggest issue is that functions can return
 * NULL, and we need to present that to the user in some reasonable
 * way that won't be confused with an ordinary return -- this is why
 * we use 0 to represent a null return for a string or an object, but
 * we use an empty array to represent a null return for a numeric
 * object.  A secondary issue has to do with strings -- we won't know the
 * size of a return string until after we have it in hand.  Both because
 * of the possibility of NULL returns and because of the indeterminate size
 * of strings, we allocate space to hold return values *with* the call,
 * rather than setting up space before the call.
 */

void mex_make_call(FILE* fp, Var* v, int first)
{
    if (!v)
        return;
    if (!first)
        fprintf(fp, ", ");

    char namebuf[128];
    if (v->tinfo == VT_obj || v->tinfo == VT_r_obj)
        fprintf(fp, "*%s", vname(v, namebuf));
    else if (v->tinfo == VT_mx && v->iospec == 'o')
        fprintf(fp, "plhs+%d", v->output_label);
    else if (v->tinfo == VT_p_scalar ||
             v->tinfo == VT_p_cscalar ||
             v->tinfo == VT_p_zscalar)
        fprintf(fp, "&%s", vname(v, namebuf));
    else if (v->tinfo == VT_const)
        fprintf(fp, "%s", v->name);
    else
        fprintf(fp, "%s", vname(v, namebuf));
    mex_make_call(fp, v->next, 0);
}


void mex_make_call(FILE* fp, Func* f)
{
    if (f->thisv)
        fprintf(fp, "in0_->");
    if (strcmp(f->funcv, "new") == 0)
        fprintf(fp, "new %s(", f->classv);
    else {
        if (f->fort)
            fprintf(fp, "MWF77_");
        fprintf(fp, "%s(", f->funcv);
    }
    mex_make_call(fp, f->args, 1);
    fprintf(fp, ")");
}


void mex_make_stmt(FILE* fp, Func* f)
{
    if (!f)
        return;

    if (f->thisv)
        fprintf(fp, 
                "    if (!in0_) {\n"
                "        mw_err_txt_ = \"Cannot dispatch to NULL\";\n"
                "        goto mw_err_label;\n"
                "    }\n");

    if (mw_generate_catch)
        fprintf(fp, "    try {\n    ");

    if (f->ret) {
        Var* v = f->ret;
        if (v->tinfo == VT_obj) {

            if (is_mxarray_type(v->basetype)) {
                fprintf(fp, "    plhs[0] = mxWrapSet_%s(&(", v->basetype);
                mex_make_call(fp, f);
                fprintf(fp ,"));\n");
            } else {
                fprintf(fp, "    out0_ = new %s(", v->basetype);
                mex_make_call(fp, f);
                fprintf(fp ,");\n");
            }

        } else if (is_array(v->tinfo)) {

            fprintf(fp, "    plhs[0] = mxWrapReturn_%s(", v->basetype);
            mex_make_call(fp, f);
            fprintf(fp, ", ");
            Expr* e = v->qual->args;
            if (e && e->next && !e->next->next) {
                fprintf(fp, " dim%d_, dim%d_);\n", 
                        e->input_label, e->next->input_label);
            } else {
                mex_alloc_size(fp, v->qual->args);
                fprintf(fp, ", 1);\n");
            }

        } else if (v->tinfo == VT_scalar  || v->tinfo == VT_r_scalar  ||
                   v->tinfo == VT_cscalar || v->tinfo == VT_r_cscalar ||
                   v->tinfo == VT_zscalar || v->tinfo == VT_r_zscalar) {

            fprintf(fp, "    out0_ = ");
            mex_make_call(fp, f);
            fprintf(fp, ";\n");

        } else if (v->tinfo == VT_string) {

            fprintf(fp, "    plhs[0] = mxWrapStrncpy(");
            mex_make_call(fp, f);
            fprintf(fp, ");\n");

        } else if (v->tinfo == VT_mx) {

            fprintf(fp, "    plhs[0] = ");
            mex_make_call(fp, f);
            fprintf(fp, ";\n");

        } else if (v->tinfo == VT_p_obj) {

            if (is_mxarray_type(v->basetype)) {
                fprintf(fp, "    plhs[0] = mxWrapSet_%s(", v->basetype);
                mex_make_call(fp, f);
                fprintf(fp, ");\n");
            } else {
                fprintf(fp, "    out0_ = ");
                mex_make_call(fp, f);
                fprintf(fp, ";\n");
            }

        } else if (v->tinfo == VT_p_scalar ||
                   v->tinfo == VT_p_cscalar ||
                   v->tinfo == VT_p_zscalar) {

            fprintf(fp, "    plhs[0] = mxWrapReturn_%s(", v->basetype);
            mex_make_call(fp, f);
            fprintf(fp, ", 1, 1);\n");

        } else if (v->tinfo == VT_r_obj) {

            if (is_mxarray_type(v->basetype)) {
                fprintf(fp, "    plhs[0] = mxWrapSet_%s(&(", v->basetype);
                mex_make_call(fp, f);
                fprintf(fp, "));\n");
            } else {
                fprintf(fp, "    out0_ = &(");
                mex_make_call(fp, f);
                fprintf(fp, ");\n");
            }

        }
    } else {
        fprintf(fp, "    ");
        mex_make_call(fp, f);
        fprintf(fp, ";\n");
    }

    if (mw_generate_catch)
        fprintf(fp, 
                "    } catch(...) {\n"
                "        mw_err_txt_ = \"Caught C++ exception from %s\";\n"
                "    }\n"
                "    if (mw_err_txt_)\n"
                "        goto mw_err_label;\n", 
                f->funcv);
}


/* -- Marshal the results -- */
/*
 * Copy inout arguments, output arguments, and any simple return
 * arguments into the plhs array.
 */


void mex_marshal_array(FILE* fp, Var* v)
{
    Expr* e = v->qual->args;
    char namebuf[128];
    const char* mtype = complex_tinfo(v) ? "mxCOMPLEX" : "mxREAL";
    const char* ws = "    ";

    if (v->tinfo == VT_rarray) {
        ws = "        ";
        fprintf(fp, "    if (out%d_ == NULL) {\n", v->output_label);
        fprintf(fp, "        plhs[%d] = mxCreateDoubleMatrix(0,0, mxREAL);\n",
                v->output_label);
        fprintf(fp, "    } else {\n");
    }

    if (!e) {

        // No dimension info -- must be an inout array
        fprintf(fp, "%splhs[%d] = mxCreateDoubleMatrix("
                "mxGetM(prhs[%d]), mxGetN(prhs[%d]), %s);\n", ws, 
                v->output_label, v->input_label, v->input_label, mtype);
        fprintf(fp, "%smxWrapCopy_%s(plhs[%d], in%d_, ", ws, 
                v->basetype, v->output_label, v->input_label);
        fprintf(fp, "mxGetM(prhs[%d])*mxGetN(prhs[%d])",
                v->input_label, v->input_label);
        fprintf(fp, ");\n");

    } else if (!(e->next)) {

        // Only size given
        fprintf(fp, "%splhs[%d] = mxCreateDoubleMatrix(dim%d_, 1, %s);\n", ws,
                v->output_label, e->input_label, mtype);
        fprintf(fp, "%smxWrapCopy_%s(plhs[%d], %s, ", ws, 
                v->basetype, v->output_label, vname(v, namebuf));
        fprintf(fp, "dim%d_",
                e->input_label);
        fprintf(fp, ");\n");

    } else if (!(e->next->next)) {

        // Two dimensions given
        fprintf(fp, "%splhs[%d] = mxCreateDoubleMatrix("
                "dim%d_, dim%d_, %s);\n", ws, 
                v->output_label, e->input_label, e->next->input_label, mtype);
        fprintf(fp, "%smxWrapCopy_%s(plhs[%d], %s, ", ws, 
                v->basetype, v->output_label, vname(v, namebuf));
        fprintf(fp, "dim%d_*dim%d_",
                e->input_label, e->next->input_label);
        fprintf(fp, ");\n");

    } else {

        // Three or more dimensions given -- punt and make it 1D
        fprintf(fp, "%splhs[%d] = mxCreateDoubleMatrix(", ws, 
                v->output_label);
        mex_alloc_size(fp, e);
        fprintf(fp, ", 1, %s);\n", mtype);
        fprintf(fp, "%smxWrapCopy_%s(plhs[%d], %s, ", ws, 
                v->basetype, v->output_label, vname(v, namebuf));
        mex_alloc_size(fp, e);
        fprintf(fp, ");\n");

    }

    if (v->tinfo == VT_rarray) {
        fprintf(fp, "    }\n");
    }

}


void mex_marshal_result(FILE* fp, Var* v, bool return_flag)
{
    char namebuf[128];
    if (is_obj(v->tinfo) && is_mxarray_type(v->basetype)) {
        if (!return_flag)
            fprintf(fp, "    plhs[%d] = mxWrapSet_%s(%s);\n",
                    v->output_label, v->basetype, vname(v,namebuf));
    } else if (is_obj(v->tinfo))
        fprintf(fp, "    plhs[%d] = mxWrapCreateP(out%d_, \"%s:%%p\");\n",
                v->output_label, v->output_label, v->basetype);
    else if (is_array(v->tinfo) ||
             v->tinfo == VT_rarray)
        mex_marshal_array(fp, v);
    else if (v->tinfo == VT_scalar ||
             v->tinfo == VT_r_scalar ||
             v->tinfo == VT_p_scalar)
        fprintf(fp, 
                "    plhs[%d] = mxCreateDoubleMatrix(1, 1, mxREAL);\n"
                "    *mxGetPr(plhs[%d]) = %s;\n",
                v->output_label, v->output_label, vname(v, namebuf));
    else if (v->tinfo == VT_cscalar || v->tinfo == VT_zscalar || 
             v->tinfo == VT_r_cscalar || v->tinfo == VT_r_zscalar || 
             v->tinfo == VT_p_cscalar || v->tinfo == VT_p_zscalar)
        fprintf(fp, 
                "    plhs[%d] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);\n"
                "    *mxGetPr(plhs[%d]) = real_%s(%s);\n"
                "    *mxGetPi(plhs[%d]) = imag_%s(%s);\n",
                v->output_label, 
                v->output_label, v->basetype, vname(v, namebuf),
                v->output_label, v->basetype, vname(v, namebuf));
    else if (v->tinfo == VT_string)
        fprintf(fp, "    plhs[%d] = mxCreateString(%s);\n",
                v->output_label, vname(v, namebuf));
}


void mex_marshal_results(FILE* fp, Var* v, bool return_flag)
{
    if (!v)
        return;
    if (v->iospec != 'i')
        mex_marshal_result(fp, v, return_flag);
    mex_marshal_results(fp, v->next, return_flag);
}


void mex_marshal_results(FILE* fp, Func* f)
{
    if (!nullable_return(f))
        mex_marshal_results(fp, f->ret, true);
    mex_marshal_results(fp, f->args, false);
}


/* -- Deallocate temporary data -- */
/*
 * Free scratch arrays allocated to hold input copies and output
 * storage.
 */


void mex_dealloc(FILE* fp, Var* v, bool return_flag)
{
    if (!v)
        return;

    if (is_array(v->tinfo) || v->tinfo == VT_string) {
        if (v->iospec == 'o')
            fprintf(fp, "    if (out%d_) mxFree(out%d_);\n", 
                    v->output_label, v->output_label);
        else if (v->iospec == 'b' || strcmp(v->basetype, "double") != 0)
            fprintf(fp, "    if (in%d_)  mxFree(in%d_);\n", 
                    v->input_label, v->input_label);
    } else if (is_obj(v->tinfo) && is_mxarray_type(v->basetype)) {
        if (v->iospec == 'i' || v->iospec == 'b')
            fprintf(fp, "    if (in%d_)  mxWrapFree_%s(in%d_);\n",
                    v->input_label, v->basetype, v->input_label);
        else if (v->iospec == 'o' && !return_flag)
            fprintf(fp, "    if (out%d_) mxWrapFree_%s(out%d_);\n",
                    v->output_label, v->basetype, v->output_label);
    }

    mex_dealloc(fp, v->next, return_flag);
}


void mex_dealloc(FILE* fp, Func* f)
{
    if (!nullable_return(f))
        mex_dealloc(fp, f->ret, true);
    mex_dealloc(fp, f->args, false);
}


/* -- Generate mex stub -- */
/*
 * Print a MEX stub function.  We generate one of these for each
 * distinct call line in the input files.
 */


void print_c_comment(FILE* fp, Func* f)
{
    fprintf(fp, "/* ---- %s: %d ----\n * ", f->fname.c_str(), f->line);
    print(fp, f);
    for (Func* fsame = f->same_next; fsame; fsame = fsame->next)
        fprintf(fp, " * Also at %s: %d\n", fsame->fname.c_str(), fsame->line);
    fprintf(fp, " */\n");
}


void print_mex_stub(FILE* fp, Func* f)
{
    print_c_comment(fp, f);
    fprintf(fp,
            "const char* stubids%d_ = \"%s\";\n\n"
            "void mexStub%d(int nlhs, mxArray* plhs[],\n"
            "              int nrhs, const mxArray* prhs[])\n"
            "{\n"
            "    const char* mw_err_txt_ = 0;\n", 
            f->id, id_string(f).c_str(), f->id);

    mex_declare_args(fp, f);
    mex_unpack_dims(fp, f);
    mex_check_dims(fp, f);
    mex_unpack_inputs(fp, f);
    mex_check_inputs(fp, f);
    mex_alloc_output(fp, f);
    mex_record_call(fp, f);
    mex_make_stmt(fp, f);
    mex_marshal_results(fp, f);
    fprintf(fp, "\nmw_err_label:\n");
    mex_dealloc(fp, f);
    fprintf(fp, 
            "    if (mw_err_txt_)\n"
            "        mexErrMsgTxt(mw_err_txt_);\n"
            "}\n\n");
}


/* -- Generate profiler output routine -- */
/*
 * The profiler / coverage testing routines include starting,
 * ending, and reporting.  Starting the profiler locks the MEX file
 * in memory.
 */

void make_profile_output(FILE* fp, Func* f, const char* printfunc)
{
    fprintf(fp,
            "        if (!mexprofrecord_)\n"
            "            %s\"Profiler inactive\\n\");\n",
            printfunc);
    for (; f; f = f->next) {
        fprintf(fp, "        %s\"%%d calls to %s:%d", 
                printfunc, f->fname.c_str(), f->line);
        for (Func* fsame = f->same_next; fsame; fsame = fsame->next)
            fprintf(fp, " (%s:%d)", fsame->fname.c_str(), fsame->line);
        fprintf(fp, "\\n\", mexprofrecord_[%d]);\n", f->id);
    }
}


/* -- Generate MEX file -- */
/*
 * Generate the overall C file to be fed into MEX.  This consists of
 * basic support code (from mwrap-support.c), type conversion
 * routines, all the call stubs, and a main dispatcher that decides
 * which stub to call.
 */


void print_mex_stubs(FILE* fp, Func* f)
{
    for (; f; f = f->next)
        print_mex_stub(fp, f);
}


void print_mex_else_cases(FILE* fp, Func* f)
{
    for (Func* fcall = f; fcall; fcall = fcall->next)
        fprintf(fp, 
                "    else if (strcmp(id, stubids%d_) == 0)\n"
                "        mexStub%d(nlhs,plhs, nrhs-1,prhs+1);\n",
                fcall->id, fcall->id);

    int maxid = max_routine_id(f);
    fprintf(fp, 
            "    else if (strcmp(id, \"*profile on*\") == 0) {\n"
            "        if (!mexprofrecord_) {\n"
            "            mexprofrecord_ = (int*) malloc(%d * sizeof(int));\n"
            "            mexLock();\n"
            "        }\n"
            "        memset(mexprofrecord_, 0, %d * sizeof(int));\n"
            "    } else if (strcmp(id, \"*profile off*\") == 0) {\n"
            "        if (mexprofrecord_) {\n"
            "            free(mexprofrecord_);\n"
            "            mexUnlock();\n"
            "        }\n"
            "        mexprofrecord_ = NULL;\n"
            "    } else if (strcmp(id, \"*profile report*\") == 0) {\n",
            maxid+1, maxid+1);

    make_profile_output(fp, f, "mexPrintf(");

    fprintf(fp,
            "    } else if (strcmp(id, \"*profile log*\") == 0) {\n"
            "        FILE* logfp;\n"
            "        if (nrhs != 2 || mxGetString(prhs[1], id, sizeof(id)) != 0)\n"
            "            mexErrMsgTxt(\"Must have two string arguments\");\n"
            "        logfp = fopen(id, \"w+\");\n"
            "        if (!logfp)\n"
            "            mexErrMsgTxt(\"Cannot open log for output\");\n");

    make_profile_output(fp, f, "fprintf(logfp, ");

    fprintf(fp, "        fclose(logfp);\n");

    fprintf(fp, 
            "    } else\n"
            "        mexErrMsgTxt(\"Unknown identifier\");\n");
}


const char* mwrap_banner = 
    "/* --------------------------------------------------- */\n"
    "/* Automatically generated by mwrap                    */\n"
    "/* --------------------------------------------------- */\n\n";

const char* mexBase =
    "/* ----\n"
    " */\n"
    "void mexFunction(int nlhs, mxArray* plhs[],\n"
    "                 int nrhs, const mxArray* prhs[])\n"
    "{\n"
    "    char id[512];\n"
    "    if (nrhs == 0) {\n"
    "        mexPrintf(\"Mex function installed\\n\");\n"
    "        return;\n"
    "    }\n\n"
    "    if (mxGetString(prhs[0], id, sizeof(id)) != 0)\n"
    "        mexErrMsgTxt(\"Identifier should be a string\");\n";

void print_mex_file(FILE* fp, Func* f)
{
    fprintf(fp, "%s", mwrap_banner);
    fprintf(fp, "%s", mex_header);

    if (mw_use_c99_complex)
        mex_c99_complex(fp);
    else if (mw_use_cpp_complex)
        mex_cpp_complex(fp);

    mex_define_copiers(fp);
    mex_casting_getters(fp);

    if (has_fortran(f)) {
        mex_define_fnames(fp, f);
        mex_fortran_decls(fp, f);
    }

    print_mex_stubs(fp, f);
    fprintf(fp, "%s", mexBase);
    print_mex_else_cases(fp, f);
    fprintf(fp, "}\n\n");
}
