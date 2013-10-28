%{
/*
 * mwrap.y
 *   Parser for mwrap.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include <stdlib.h>
#include <string.h>
#include <string>
#include "mwrap-ast.h"

extern "C" {
    int yylex();
    int yywrap();
    int yyerror(const char* s);
}

using std::string;

bool  mw_generate_catch = false;  // Catch C++ exceptions?
bool  mw_use_cpp_complex = false; // Use C++ complex types?
bool  mw_use_c99_complex = false; // Use C99 complex types?
bool  mw_promote_int = false;     // Convert integer types to mwSize?
int   listing_flag = 0;           // Output filenames from @ commands?
int   mbatching_flag = 0;         // Output on @ commands?
int   linenum = 0;                // Lexer line number
FILE* outfp   = 0;                // MATLAB output file
FILE* outcfp  = 0;                // C output file

static int    type_errs = 0;            // Number of typecheck errors
static int    func_id = 0;              // Assign stub numbers
static Func*  funcs   = 0;              // AST - linked list of functions
static Func*  lastfunc = 0;             // Last link in funcs list
static const char*  mexfunc = "mexfunction";  // Name of mex function
static string current_ifname;           // Current input file name


#define MAX_INCLUDE_DEPTH 10
static string include_stack_names[MAX_INCLUDE_DEPTH];
extern int include_stack_ptr;

extern "C" void set_include_name(const char* s)
{
    include_stack_names[include_stack_ptr] = current_ifname;
    current_ifname = s;
}

extern "C" void get_include_name()
{
    current_ifname = include_stack_names[include_stack_ptr].c_str();
}


inline void add_func(Func* func)
{
    static std::map<string,Func*> func_lookup;
    if (!funcs) {
        funcs = func;
        lastfunc = func;
        return;
    } 

    Func*& func_ptr = func_lookup[id_string(func)];
    if (func_ptr) {
        func_ptr->same_next = func;
    } else {
        lastfunc->next = func;
        lastfunc = func;
    }
    func_ptr = func;
}

%}

%union {
    char* string;
    struct Func* func;
    struct Var* var;
    struct TypeQual* qual;
    struct Expr* expr;
    struct InheritsDecl* inherits;
    char c;
}

%token NON_C_LINE
%token NEW TYPEDEF CLASS FORTRAN
%token <string> ID 
%token <string> NUMBER STRING
%token <char> INPUT OUTPUT INOUT

%type <func> func funcall
%type <var>  var basevar args argsrest
%type <c>    iospec
%type <qual> quals aqual
%type <expr> arrayspec exprs exprrest expr
%type <inherits> inheritslist inheritsrest

%error-verbose

%%
statements: statement statements | ;

statement:
  basevar '=' funcall { 
      $3->ret = $1; 
      $3->id = ++func_id;
      type_errs += typecheck($3, linenum);
      if (outfp)
          print_matlab_call(outfp, $3, mexfunc); 
      add_func($3);
  }
  | funcall { 
      $1->id = ++func_id;
      type_errs += typecheck($1, linenum);
      if (outfp)
          print_matlab_call(outfp, $1, mexfunc); 
      add_func($1);
  } 
  | tdef 
  | classdef 
  | NON_C_LINE 
  | error ';' { yyerrok; } ;

tdef: 
  TYPEDEF ID ID ';' { 
      if (strcmp($2, "numeric") == 0) {
          add_scalar_type($3);
      } else if (strcmp($2, "dcomplex") == 0) {
          add_zscalar_type($3);
      } else if (strcmp($2, "fcomplex") == 0) {
          add_cscalar_type($3);
      } else if (strcmp($2, "mxArray") == 0) {
          add_mxarray_type($3);
      } else {
          fprintf(stderr, "Unrecognized typespace: %s\n", $2);
          ++type_errs;
      }
      delete[] $2;
      delete[] $3;
  } ;

classdef:
  CLASS ID ':' inheritslist ';' {
      add_inherits($2, $4);
      delete[] $2;
      destroy($4);
  }

inheritslist:
  ID inheritsrest { $$ = new InheritsDecl($1, $2); } ;

inheritsrest:
  ',' ID inheritsrest { $$ = new InheritsDecl($2, $3); }
  | { $$ = NULL; } ;

funcall: func '(' args ')' ';' { $$ = $1; $$->args = $3; } ;

args:
  var argsrest { $$ = $1; $$->next = $2; }
  | { $$ = NULL; } ;

argsrest:
  ',' var argsrest {$$ = $2; $$->next = $3; }
  | { $$ = NULL; } ;

basevar: ID ID               { $$ = new Var('o', promote_int($1), NULL, $2); }
basevar: ID quals ID         { $$ = new Var('o', promote_int($1), $2,   $3); }
basevar: ID ID aqual         { $$ = new Var('o', promote_int($1), $3,   $2); }

var: iospec ID ID            { $$ = new Var($1,  promote_int($2), NULL, $3); }
var: iospec ID quals ID      { $$ = new Var($1,  promote_int($2), $3,   $4); }
var: iospec ID ID aqual      { $$ = new Var($1,  promote_int($2), $4,   $3); }

var: iospec ID NUMBER        { $$ = new Var($1,  promote_int($2), NULL, $3); }
var: iospec ID quals NUMBER  { $$ = new Var($1,  promote_int($2), $3,   $4); }

var: iospec ID STRING        { $$ = new Var($1,  promote_int($2), NULL, $3); }
var: iospec ID quals STRING  { $$ = new Var($1,  promote_int($2), $3,   $4); }

iospec: 
  INPUT    { $$ = 'i'; }
  | OUTPUT { $$ = 'o'; }
  | INOUT  { $$ = 'b'; }
  |        { $$ = 'i'; } ;

quals: 
  '*'         { $$ = new TypeQual('*', NULL); }
  | '&'       { $$ = new TypeQual('&', NULL); } 
  | aqual     { $$ = $1; } ;

aqual:
  arrayspec       { $$ = new TypeQual('a', $1); } 
  | arrayspec '&' { $$ = new TypeQual('r', $1); } ;

arrayspec: '[' exprs ']' { $$ = $2; } ;

exprs: 
  expr exprrest { $$ = $1; $$->next = $2; }
  |             { $$ = NULL; } 

exprrest: 
  ',' expr exprrest { $$ = $2; $$->next = $3; }
  |                 { $$ = NULL; }

expr: 
  ID       { $$ = new Expr($1); }
  | NUMBER { $$ = new Expr($1); }

func: 
  ID '-' '>' ID '.' ID { $$ = new Func($1, $4, $6, current_ifname, linenum); }
  | ID          { $$ = new Func(NULL, NULL, $1, current_ifname, linenum); } 
  | FORTRAN ID  { $$ = new Func(NULL, NULL, $2, current_ifname, linenum); 
                  $$->fort = true;
                } 
  | NEW ID  { $$ = new Func(NULL, $2, mwrap_strdup("new"), 
                          current_ifname, linenum); 
            }
  ;

%%
#include <stdio.h>
#include <string.h>

extern FILE* yyin;

int yywrap()
{
    return 1;
}

int yyerror(const char* s)
{
    fprintf(stderr, "Parse error (%s:%d): %s\n", current_ifname.c_str(),
            linenum, s);
}

char* mwrap_strdup(const char* s)
{
    char* result = new char[strlen(s)+1];
    strcpy(result, s);
    return result;
}

const char* help_string = 
"mwrap 0.33.3 - MEX file generator for MATLAB and Octave\n"
"\n"
"Syntax:\n"
"  mwrap [-mex outputmex] [-m output.m] [-c outputmex.c] [-mb]\n"
"        [-list] [-catch] infile1 infile2 ...\n"
"\n"
"  -mex outputmex -- specify the MATLAB mex function name\n"
"  -m output.m    -- generate the MATLAB stub called output.m\n"
"  -c outputmex.c -- generate the C file outputmex.c\n"
"  -mb            -- generate .m files specified with @ redirections\n"
"  -list          -- list files specified with @ redirections\n"
"  -catch         -- generate C++ exception handling code\n"
"  -im            -- convert int, long, uint, and ulong types to mwSize\n"
"  -c99complex    -- add support code for C99 complex types\n"
"  -cppcomplex    -- add support code for C++ complex types\n"
"\n";

int main(int argc, char** argv)
{
    int j;
    int err_flag = 0;
    init_scalar_types();

    if (argc == 1) {
        fprintf(stderr, "%s", help_string);
        return 0;
    } else {
        for (j = 1; j < argc; ++j) {
            if (strcmp(argv[j], "-m") == 0 && j+1 < argc)
                outfp = fopen(argv[j+1], "w+");
            if (strcmp(argv[j], "-c") == 0 && j+1 < argc)
                outcfp = fopen(argv[j+1], "w+");
            if (strcmp(argv[j], "-mex") == 0 && j+1 < argc)
                mexfunc = argv[j+1];
            if (strcmp(argv[j], "-mb") == 0)
                mbatching_flag = 1;
            if (strcmp(argv[j], "-list") == 0)
                listing_flag = 1;
            if (strcmp(argv[j], "-catch") == 0)
                mw_generate_catch = true;
            if (strcmp(argv[j], "-im") == 0)
                mw_promote_int = true;
            if (strcmp(argv[j], "-c99complex") == 0) 
                mw_use_c99_complex = true;
            if (strcmp(argv[j], "-cppcomplex") == 0) 
                mw_use_cpp_complex = true;
        }

        if (mw_use_c99_complex || mw_use_cpp_complex) {
            add_zscalar_type("dcomplex");
            add_cscalar_type("fcomplex");
        }

        for (j = 1; j < argc; ++j) {
            if (strcmp(argv[j], "-m") == 0 ||
                strcmp(argv[j], "-c") == 0 ||
                strcmp(argv[j], "-mex") == 0)
                ++j;
            else if (strcmp(argv[j], "-mb") == 0 ||
                     strcmp(argv[j], "-list") == 0 ||
                     strcmp(argv[j], "-catch") == 0 ||
                     strcmp(argv[j], "-im") == 0 ||
                     strcmp(argv[j], "-c99complex") == 0 ||
                     strcmp(argv[j], "-cppcomplex") == 0);
            else {
                linenum = 1;
                type_errs = 0;
                yyin = fopen(argv[j], "r");
                if (yyin) {
                    current_ifname = argv[j];
                    err_flag += yyparse();
                    fclose(yyin);
                } else {
                    fprintf(stderr, "Could not read %s\n", argv[j]);
                }
                if (type_errs)
                    fprintf(stderr, "%s: %d type errors detected\n", 
                            argv[j], type_errs);
                err_flag += type_errs;
            }
        }
    }
    if (!err_flag && outcfp)
        print_mex_file(outcfp, funcs);
    destroy(funcs);
    destroy_inherits();
    if (outfp)
        fclose(outfp);
    if (outcfp)
        fclose(outcfp);
    return err_flag;
}
