/*
 * mwrap-ast.cc
 *   Recursive routines for AST printing, identifier construction,
 *   and destruction.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cassert>
#include "mwrap-ast.h"


/* -- Add inheritance relationship -- */
/*
 * We represent inheritance relationships as a list of child classes
 * associated with each parent.  This makes sense when we generate the
 * type casting routines, since we only allow casts to the parent
 * types.  But it's the opposite of how inheritance relations are
 * specified syntactically.
 */

map<string, InheritsDecl*> class_decls;


void add_inherits(const char* childname, InheritsDecl* ilist)
{
    for (; ilist; ilist = ilist->next) {
        InheritsDecl*& children = class_decls[ilist->name];
        children = new InheritsDecl(mwrap_strdup(childname), children);
    }
}


/* -- Add scalar type data -- */
/*
 * The names of the numeric scalar types are stored in the set scalar_decls.
 * We provide a basic list, but the user can specify more -- the only real
 * constraint is there has to be a meaningful typecast to/from a double.
 */

set<string> scalar_decls;
set<string> cscalar_decls;
set<string> zscalar_decls;
set<string> mxarray_decls;


void init_scalar_types()
{
    const char* scalar_types[] = {
        "double", "float", 
        "long", "int", "char", 
        "ulong", "uint", "uchar",
        "bool", "size_t", 
        "mwSize", "mwIndex", "mwSignedIndex", NULL};

    for (const char** s = scalar_types; *s; ++s)
        add_scalar_type(*s);
}


char *promote_int(char* name)
{
  if( mw_promote_int == true )
    {      
      if( strcmp(name,"int") == 0 ) return strdup("mwSize");
      if( strcmp(name,"long") == 0 ) return strdup("mwSize");
      if( strcmp(name,"uint") == 0 ) return strdup("mwSize");
      if( strcmp(name,"ulong") == 0 ) return strdup("mwSize");
    }
  return name;
}



void add_scalar_type(const char* name)
{
    scalar_decls.insert(name);
}


void add_cscalar_type(const char* name)
{
    cscalar_decls.insert(name);
}


void add_zscalar_type(const char* name)
{
    zscalar_decls.insert(name);
}


void add_mxarray_type(const char* name)
{
    mxarray_decls.insert(name);
}


bool is_scalar_type(const char* name)
{
    return (scalar_decls.find(name) != scalar_decls.end());
}


bool is_cscalar_type(const char* name)
{
    return (cscalar_decls.find(name) != cscalar_decls.end());
}


bool is_zscalar_type(const char* name)
{
    return (zscalar_decls.find(name) != zscalar_decls.end());
}


bool is_mxarray_type(const char* name)
{
    return (mxarray_decls.find(name) != mxarray_decls.end());
}


/* -- Print the AST -- */
/*
 * Print a human-readable translation of the abstract syntax tree to
 * the output.  This is only used for generating comprehensible comments
 * in the C code.
 */


void print(FILE* fp, Expr* e)
{
    if (!e)
        return;
    fprintf(fp, "%s", e->value);
    if (e->next) {
        fprintf(fp, ", ");
        print(fp, e->next);
    }
}


void print(FILE* fp, TypeQual* q)
{
    if (!q)
        return;
    if (q->qual == 'a') {
        fprintf(fp, "[");
        print(fp, q->args);
        fprintf(fp, "]");
    } else
        fprintf(fp, "%c", q->qual);
}


void print_iospec(FILE* fp, Var* v)
{
    if (v->iospec == 'o')
        fprintf(fp, "output ");
    else if (v->iospec == 'b')
        fprintf(fp, "inout ");
}


void print_var(FILE* fp, Var* v)
{
    fprintf(fp, "%s", v->basetype);
    print(fp, v->qual);
    fprintf(fp, " %s", v->name);
}


void print_args(FILE* fp, Var* v)
{
    if (!v)
        return;
    print_iospec(fp, v);
    print_var(fp, v);
    if (v->next) {
        fprintf(fp, ", ");
        print_args(fp, v->next);
    }
}


void print(FILE* fp, Func* f)
{
    if (!f)
        return;
    if (f->ret) {
        print_var(fp, f->ret);
        fprintf(fp, " = ");
    }
    if (f->thisv)
        fprintf(fp, "%s->%s.", f->thisv, f->classv);
    fprintf(fp, "%s(", f->funcv);
    print_args(fp, f->args);
    fprintf(fp, ");\n");
}


/* -- ID string -- */
/*
 * The ID string is a compressed text representation of the AST for
 * a function.  We replace all the variable names by the letter 'x'
 * so that when we have the same function call with different arguments,
 * we will generate the same signature.
 */


string id_string(Expr* e)
{
    if (!e)
        return "";
    return "x" + id_string(e->next);
}


string id_string(TypeQual* q)
{
    if (!q)
        return "";
    if (q->qual == 'a')
        return "[" + id_string(q->args) + "]";
    else
        return (string() + q->qual);
}


string id_string(Var* v)
{
    if (!v)
        return "";

    string name;
    if (v->iospec == 'i')
        name += "i ";
    else if (v->iospec == 'o')
        name += "o ";
    else
        name += "io ";
    name += v->basetype;
    name += id_string(v->qual);
    if (v->tinfo == VT_const) {
        name += " ";
        name += v->name;
    }
    if (v->next) {
        name += ", ";
        name += id_string(v->next);
    }
    return name;
}


string id_string(Func* f)
{
    if (!f)
        return "";

    string name;
    if (f->ret) {
        name += id_string(f->ret);
        name += " = ";
    }
    if (f->thisv) {
        name += f->thisv;
        name += "->";
        name += f->classv;
        name += ".";
    }
    name += f->funcv;
    name += "(";
    name += id_string(f->args);
    name += ")";
    return name;
}


/* -- Destroy an AST -- */


void destroy(Expr* e)
{
    if (!e)
        return;
    destroy(e->next);
    delete[] e->value;
    delete e;
}


void destroy(TypeQual* q)
{
    if (!q)
        return;
    destroy(q->args);
    delete q;
}


void destroy(Var* v)
{
    if (!v) 
        return;
    destroy(v->next);
    destroy(v->qual);
    delete[] v->basetype;
    delete[] v->name;
    delete v;
}


void destroy(Func* f)
{
    while (f) {
        Func* oldf = f;
        f = f->next;

        destroy(oldf->same_next);
        destroy(oldf->ret);
        if (oldf->thisv)
            delete[] oldf->thisv;
        if (oldf->classv)
            delete[] oldf->classv;
        delete[] oldf->funcv;
        destroy(oldf->args);
        delete oldf;
    }
}


void destroy(InheritsDecl* i)
{
    if (!i)
        return;
    destroy(i->next);
    delete i;
}


void destroy_inherits()
{
    map<string,InheritsDecl*>::iterator i = class_decls.begin();
    for (; i != class_decls.end(); ++i) {
        destroy(i->second);
        i->second = NULL;
    }
}
