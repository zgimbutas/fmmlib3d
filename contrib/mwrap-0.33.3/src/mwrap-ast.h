/*
 * mwrap-ast.h
 *   MWrap abstract syntax tree declarations.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#ifndef MWRAP_AST_H
#define MWRAP_AST_H

#include <stdio.h>
#include <string>
#include <map>
#include <set>

using std::string;
using std::map;
using std::set;

enum {
    VT_unk,         // Unknown
    VT_obj,         // Object
    VT_array,       // Numeric array (real)
    VT_carray,      // Numeric array (single complex)
    VT_zarray,      // Numeric array (double complex)
    VT_rarray,      // Reference to numeric array
    VT_scalar,      // Numeric scalar (real)
    VT_cscalar,     // Numeric scalar (single complex)
    VT_zscalar,     // Numeric scalar (double complex)
    VT_string,      // Char string
    VT_mx,          // mxArray
    VT_p_obj,       // Pointer to object
    VT_p_scalar,    // Pointer to scalar (real)
    VT_p_cscalar,   // Pointer to scalar (single complex)
    VT_p_zscalar,   // Pointer to scalar (double complex)
    VT_r_obj,       // Reference to object
    VT_r_scalar,    // Reference to scalar (real)
    VT_r_cscalar,   // Reference to scalar (single complex)
    VT_r_zscalar,   // Reference to scalar (double complex)
    VT_const,       // Constant expression
};

struct InheritsDecl {
    InheritsDecl(char* name, InheritsDecl* next) : 
        name(name), next(next) {}

    char* name;              // Name of inherited class
    InheritsDecl* next;
};

struct Expr {
    Expr(char* value) : value(value), next(NULL), input_label(-1) {}

    int input_label;  // Index of dim variable in input arg list
    char* value;      // Expression to be evaluated to get dimension
    Expr* next;
};

struct TypeQual {
    TypeQual(char qual, Expr* args) : 
        qual(qual), args(args) {}

    char qual;        // Pointer ('*'), ref ('&'), array ('a'), or nothing (0)
    Expr* args;       // Array / cstring dimension list
};

struct Var {
    Var(char iospec, char* basetype, TypeQual* qual, char* name) :
        iospec(iospec), basetype(basetype), qual(qual), tinfo(VT_unk),
        name(name), next(NULL), input_label(-1), output_label(-1) {}

    int input_label;  // Index in input arg list
    int output_label; // Index in output arg list
    char iospec;      // Indicate input ('i'), output ('o'), or both ('b')
    char* basetype;   // Name of the base type
    TypeQual* qual;   // Type qualifier (pointer, ref, etc)
    int tinfo;        // General type identifier (see VT_* list above)
    char* name;       // MATLAB text for variable name (or value)
    Var* next;
};

struct Func {
    Func(char* thisv, char* classv, char* funcv, 
         const string& fname, int line) :
        fname(fname), line(line), fort(false), id(-1),
        thisv(thisv), classv(classv), funcv(funcv), 
        args(NULL), ret(NULL), next(NULL), same_next(NULL) {}

    string fname; // Name of file where this function was defined
    int    line;  // Number of the line where the function was defined
    bool   fort;  // Flag whether this is a FORTRAN function

    int id;       // Identifier (used for numbering stub functions)
    char* thisv;  // This var (only for methods)
    char* classv; // Class type (only for methods or 'new' calls)
    char* funcv;  // Function or method name
    Var* args;    // Arguments
    Var* ret;     // Return variable

    Func* next;       // Next in ordinary linked list
    Func* same_next;  // Next in list of type-identical calls
};

extern "C" char* mwrap_strdup(const char* s);

extern map<string, InheritsDecl*> class_decls;
void add_inherits(const char* childname, InheritsDecl* ilist);

extern set<string> scalar_decls;
extern set<string> cscalar_decls;
extern set<string> zscalar_decls;
extern set<string> mxarray_decls;
void init_scalar_types();
char *promote_int(char* name);
void add_scalar_type(const char* name);
void add_cscalar_type(const char* name);
void add_zscalar_type(const char* name);
void add_mxarray_type(const char* name);
bool is_scalar_type(const char* name);
bool is_cscalar_type(const char* name);
bool is_zscalar_type(const char* name);
bool is_mxarray_type(const char* name);

string id_string(Func* f);
int typecheck(Func* func, int line);

void print(FILE* fp, Func* func);
void print_matlab_call(FILE* fp, Func* func, const char* mexfunc);
void print_mex_file(FILE* fp, Func* f);

void destroy(Func* func);
void destroy(InheritsDecl* ilist);
void destroy_inherits();

extern bool mw_generate_catch;
extern bool mw_use_c99_complex;
extern bool mw_use_cpp_complex;
extern bool mw_promote_int;

#endif /* MWRAP_AST_H */

