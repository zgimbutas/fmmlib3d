%{
/*
 * mwrap.l
 *   Lexer for MWrap.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include "mwrap.hh"
#include <string.h>
#include <ctype.h>

extern int listing_flag;          // Output filenames from @ commands?
extern int mbatching_flag;        // Do we want to output on @ commands?
extern int linenum;               // Lexer line number
extern FILE* outfp;               // MATLAB output file
extern FILE* outcfp;              // C output file

static int done_at_switch;        // Set when @ redirection is done

extern char* mwrap_strdup(char* s);

static int is_name_char(char c)
{
    return (isalnum(c) || c == '_');
}

static char* fname_scan_line(char* s)
{
    static char namebuf[256];  /* FIXME */
    int name_start, name_end, i, j;

    /* Name ends at last alphanum before args */
    name_end = 0;
    while (s[name_end] && s[name_end] != '(')
        ++name_end;
    while (name_end > 0 && !is_name_char(s[name_end]))
        --name_end;

    /* Back up to the start of the name */
    name_start = name_end;
    while (s[name_start] > 0 && is_name_char(s[name_start]))
        --name_start;

    /* Copy the name into the buf and add .m */
    for (i = name_start+1, j = 0; i <= name_end; ++i, ++j)
        namebuf[j] = s[i];
    strcpy(namebuf+j, ".m");
    return namebuf;
}


#define MAX_INCLUDE_DEPTH 10
YY_BUFFER_STATE include_stack[MAX_INCLUDE_DEPTH];
int include_stack_line[MAX_INCLUDE_DEPTH];
int include_stack_ptr = 0;
extern void set_include_name(const char* s);
extern void get_include_name();


/* The lexer switches states when it sees a specially formatted comment
 * as the first non-blank mode in a line.  In total, there are six states:
 *
 *   INITIAL - start of line
 *   TS      - ordinary text line
 *   AS      - at line (redirection)
 *   FS      - @function line (function declaration + redirection)
 *   SS      - embedded C line
 *   BS      - block of embedded C
 *   CS      - C call line
 *   COMMS   - Comment line
 * 
 * Transitions are as follows.
 *
 *   "$[" :        INITIAL -> BS
 *   "$]" :        BS -> INITIAL
 *   "@" :         INITIAL -> AS  (batching mode only)
 *   "@function" : INITIAL -> FS
 *   "$" :         INITIAL -> SS
 *   "#" :         INITIAL -> CS
 *   non-blank :   INITIAL -> TS
 *   newline:      (AS, SS, CS, TS) -> INITIAL
 */

%}

%s CSTATE
%s SSTATE
%s BSTATE
%s ASTATE
%s FSTATE
%s TSTATE
%s COMMSTATE
%x INCLSTATE

%%
<INITIAL>"@function" { 
    if (mbatching_flag && outfp)
        fclose(outfp);
    BEGIN FSTATE;
    return NON_C_LINE;
}

<INITIAL>"@include" { BEGIN INCLSTATE; return NON_C_LINE; }

<INITIAL>"@"   { 
    if (mbatching_flag && outfp)
        fclose(outfp); 
    done_at_switch = 0; 
    BEGIN ASTATE;
    return NON_C_LINE; 
}

<INITIAL>"#"             { BEGIN CSTATE; }
<INITIAL>\$\[[ \t\r]*\n  { BEGIN BSTATE; ++linenum; return NON_C_LINE; }
<INITIAL>"$"             { BEGIN SSTATE;            return NON_C_LINE; }
<INITIAL>"//"            { BEGIN COMMSTATE;         return NON_C_LINE; }

<INITIAL>\n    { if (outfp) fprintf(outfp, "%s", yytext); 
                 ++linenum; 
                 BEGIN 0; 
                 return NON_C_LINE; }
<INITIAL>[ \t] { if (outfp) fprintf(outfp, "%s", yytext); }
<INITIAL>.     { if (outfp) fprintf(outfp, "%s", yytext); 
                 BEGIN TSTATE; 
                 return NON_C_LINE; }

<COMMSTATE>\n  { ++linenum; BEGIN 0; }
<COMMSTATE>. ;

<<EOF>> {
    if (--include_stack_ptr < 0) {
        yyterminate();
    } else {
        yy_delete_buffer(YY_CURRENT_BUFFER);
        yy_switch_to_buffer(include_stack[include_stack_ptr]);
        linenum = include_stack_line[include_stack_ptr];
        get_include_name();
        BEGIN INCLSTATE;
    }
}

<INCLSTATE>[ \t\r]* ;
<INCLSTATE>[^ \t\r\n]+ {
    if (include_stack_ptr >= MAX_INCLUDE_DEPTH) {
        fprintf(stderr, "Error: Includes nested too deeply");
        exit(-1);
    }
    set_include_name(yytext);
    include_stack_line[include_stack_ptr] = linenum;
    include_stack[include_stack_ptr++] = YY_CURRENT_BUFFER;
    yyin = fopen(yytext, "r");
    if (!yyin) {
        fprintf(stderr, "Error: Could not read '%s'\n", yytext);
        exit(-1);
    }
    linenum = 1;
    yy_switch_to_buffer(yy_create_buffer(yyin, YY_BUF_SIZE));
    BEGIN 0;
}
<INCLSTATE>\n { ++linenum; BEGIN 0; }

<FSTATE>[^\n]+ {
    char* fname = fname_scan_line(yytext);
    if (mbatching_flag) {
        outfp = fopen(fname, "w+");
        if (!outfp) {
            fprintf(stderr, "Error: Could not write %s\n", yytext);
            exit(-1);
        }
    }
    if (listing_flag)
        fprintf(stdout, "%s\n", fname);
    if (outfp)
        fprintf(outfp, "function%s\n", yytext);
}
<FSTATE>\n { ++linenum; BEGIN 0; }

<ASTATE>[ \t\r]+ ;
<ASTATE>[^ \t\r\n]+ { if (mbatching_flag && !done_at_switch) {
                         outfp = fopen(yytext, "w+");
                         if (!outfp) {
                             fprintf(stderr, "Error: Could not write %s\n",
                                     yytext);
                             exit(-1);
                         }
                      }
                      if (listing_flag && !done_at_switch)
                          fprintf(stdout, "%s\n", yytext);
                      done_at_switch = 1; 
                    }
<ASTATE>\n          { ++linenum; BEGIN 0; }

<TSTATE>\n   { if (outfp) fprintf(outfp, "%s", yytext); ++ linenum; BEGIN 0; }
<TSTATE>.    { if (outfp) fprintf(outfp, "%s", yytext); }

<SSTATE>\n   { if (outcfp) fprintf(outcfp, "%s", yytext); ++linenum; BEGIN 0; }
<SSTATE>.    { if (outcfp) fprintf(outcfp, "%s", yytext); }

<BSTATE>\$\][ \t\r]*\n { ++linenum; BEGIN 0; }
<BSTATE>\n   { if (outcfp) fprintf(outcfp, "%s", yytext); ++linenum; }
<BSTATE>.    { if (outcfp) fprintf(outcfp, "%s", yytext); }

<CSTATE>"new"     { return NEW; }
<CSTATE>"FORTRAN" { return FORTRAN; }
<CSTATE>"input"   { return INPUT; }
<CSTATE>"output"  { return OUTPUT; }
<CSTATE>"inout"   { return INOUT; }
<CSTATE>"class"   { return CLASS; }
<CSTATE>"typedef" { return TYPEDEF; }

<CSTATE>((::)?[_a-zA-Z][_a-zA-Z0-9]*)* { 
    yylval.string = mwrap_strdup(yytext); 
    return ID; 
}
<CSTATE>[0-9]+ {
    yylval.string = mwrap_strdup(yytext); 
    return NUMBER; 
}
<CSTATE>\'[^'\n]*['\n] {
    yylval.string = mwrap_strdup(yytext);
    return STRING;
}
<CSTATE>\/\/[^\n]* ;

<CSTATE>[ \t\r]+ ;
<CSTATE>\n { ++linenum; BEGIN 0; }
<CSTATE>.   return yytext[0];

%%
