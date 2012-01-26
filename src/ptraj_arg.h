#ifndef INC_PTRAJ_ARG_H
#define INC_PTRAJ_ARG_H
#ifdef __cplusplus
extern "C" {
#endif
/*! \file ptraj_arg.h

  This file contains argStackType, the memory-safe replacement
  for the original argument stack in ptraj_actions and 
  ptraj_analyze.
 */

// ---------- Argument Stack routines ------------------------------------------
// argStackType - just used to set up arglist
typedef struct _argStackType {
  int nargs;
  char **arglist;
  char *marked;
} argStackType;

void printArgumentStack(argStackType **);
char *getArgumentString(argStackType **, char *);
int getArgumentInteger(argStackType **, int );
double getArgumentDouble(argStackType **, double );
int argumentStringContains(argStackType **, char *);
int argumentStackContains(argStackType **, char *);
int argumentStackKeyToInteger(argStackType **, char *, int );
double argumentStackKeyToDouble(argStackType **,char *, double );
float argumentStackKeyToFloat(argStackType **,char *, float );
char* argumentStackKeyToString(argStackType **, char *, char* );

#ifdef __cplusplus
}
#endif
#endif
