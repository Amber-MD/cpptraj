#ifndef INC_PTRAJMASK_H
#define INC_PTRAJMASK_H
/// PtrajMask: The enhanced atom selection mask parser from ptraj.
/// Originally written by Viktor Hornak, Stony Brook University.
/// Adapted as standalone code by Dan Roe, NIST.
// The NAME parameters control how the mask parser expects strings to look.
// Originally, the parameter file assumes the atom, residue, symbol, etc. 
// names to be four characters. When stored as a string, the NULL character
// is required, requiring a size of 5. It has been increased to 6 in cpptraj
// to accomodate the slightly larger names that can be found in MOL2 files.
#include "Name.h"
// More defines
#define  MAXSELE 1000
#define  ALL      0
#define  NUMLIST  1
#define  NAMELIST 2
#define  TYPELIST 3
#define  ELEMLIST 4
int tokenize(char *, char *);
int torpn(char *, char *);
/* parseMaskString()
 * The main interface to the mask parser. Takes a mask expression and some
 * information from a parameter file (# atoms, # residues, atom names, residue
 * names, an array containing the first atom # of each residue, atomic coords
 * in X0Y0Z0X1Y1Z1... format, atom types, and a debug value controlling how
 * much debug information is printed (internally the global int prnlev).
 * It returns a character mask array mask[i]='T'|'F', i=0,atoms-1
 * which contains the resulting atom selection
 */
char *parseMaskString(char*,int,int,NAME*,NAME*,int*,double*,NAME*,int);

#ifdef __cplusplus
extern "C"
#endif 
char* parseMaskC(char*,int,int,NAME*,NAME*,int*,double*,NAME*,int);
#endif
