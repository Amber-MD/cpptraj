#ifndef INC_PTRAJMASK_H
#define INC_PTRAJMASK_H
/*! \file PtrajMask.h
    \brief The enhanced atom selection mask parser from ptraj.
    \author Viktor Hornak, Stony Brook University.
    \author Daniel R. Roe, NIST.
 */
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
char *parseMaskString(char*,int,int,NAME*,NAME*,int*,double*,NAME*,int);

#ifdef __cplusplus
extern "C"
#endif 
char* parseMaskC(char*,int,int,NAME*,NAME*,int*,double*,NAME*,int);
#endif
