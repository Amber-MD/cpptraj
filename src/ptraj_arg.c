#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ptraj_arg.h"

// ========== ARGUMENT functions ===============================================
// Uses argStackType from ptraj_actions.h. To be consistent with original
// argument stack behavior, copies of strings are returned.

// getArgumentString()
/// Return the next unmarked argument
char *getArgumentString(argStackType **argumentStackAddress, char *defvalue) {
  int arg;
  char *retvalue = NULL;
  argStackType *argumentStack;
  argumentStack = *argumentStackAddress;
  for (arg = 0; arg < argumentStack->nargs; arg++) {
    // If arg is unmarked, mark and return copy
    if ( argumentStack->marked[arg] == 'F' ) {
      argumentStack->marked[arg] = 'T';
      retvalue = (char*) malloc( (strlen(argumentStack->arglist[arg])+1) * sizeof(char));
      strcpy(retvalue, argumentStack->arglist[arg]);
      return retvalue;
    }
  }
  // No more unmarked args; return copy of defvalue
  if (defvalue!=NULL) {
    retvalue = (char*) malloc( (strlen(defvalue)+1) * sizeof(char));
    strcpy(retvalue, defvalue);
  }
  return retvalue;
}

// getArgumentInteger()
/// Return the next unmarked argument as an integer
int getArgumentInteger(argStackType **argumentStackAddress, int defvalue) {
  int arg;
  argStackType *argumentStack;
  argumentStack = *argumentStackAddress;
  for (arg = 0; arg < argumentStack->nargs; arg++) {
    // If arg is unmarked, mark and return
    if ( argumentStack->marked[arg] == 'F' ) {
      argumentStack->marked[arg] = 'T';
      return atoi(argumentStack->arglist[arg]);
    }
  }
  // No more unmarked args; return defvalue
  return defvalue;
}

// getArgumentDouble()
/// Return the next unmarked argument as a double
double getArgumentDouble(argStackType **argumentStackAddress, double defvalue) {
  int arg;
  argStackType *argumentStack;
  argumentStack = *argumentStackAddress;
  for (arg = 0; arg < argumentStack->nargs; arg++) {
    // If arg is unmarked, mark and return
    if ( argumentStack->marked[arg] == 'F' ) {
      argumentStack->marked[arg] = 'T';
      return atof(argumentStack->arglist[arg]);
    }
  }
  // No more unmarked args; return defvalue
  return defvalue;
}

// argumentStringContains()
/// Return true and mark argument if next unmarked arg matches, false otherwise
int argumentStringContains(argStackType **argumentStackAddress, char *match) {
  int arg;
  argStackType *argumentStack;
  argumentStack = *argumentStackAddress;
  for (arg = 0; arg < argumentStack->nargs; arg++) {
    // If arg is unmarked and matches, mark and return 1, otherwise return 0
    if ( argumentStack->marked[arg] == 'F' ) {
      if (strcmp(argumentStack->arglist[arg], match)==0) {
        argumentStack->marked[arg] = 'T';
        return 1;
      } else {
        return 0;
      }
    }
  }
  return 0;
}

// argumentStackContains()
/// Return true and mark argument if argument is present, false otherwise.
int argumentStackContains(argStackType **argumentStackAddress, char *match) {
  int arg;
  argStackType *argumentStack;
  argumentStack = *argumentStackAddress;
  for (arg = 0; arg < argumentStack->nargs; arg++) {
    // If arg is unmarked and matches, mark and return 1
    if ( argumentStack->marked[arg]=='F') {
      if ( strcmp(argumentStack->arglist[arg], match)==0 ) {
        argumentStack->marked[arg] = 'T';
        return 1;
      }
    }
  }
  return 0;
}

// argumentStackKeyToInteger()
int argumentStackKeyToInteger(argStackType **argumentStackAddress, char *match, int defvalue) {
  int arg;
  argStackType *argumentStack;
  argumentStack = *argumentStackAddress;
  for (arg = 0; arg < argumentStack->nargs - 1; arg++) {
    if ( argumentStack->marked[arg]=='F' && strcmp(argumentStack->arglist[arg], match)==0 ) {
      // Check if next arg already marked
      if (argumentStack->marked[arg+1]=='T') continue; 
      argumentStack->marked[arg] = 'T';
      arg++;
      argumentStack->marked[arg] = 'T';
      return atoi(argumentStack->arglist[arg]);
    }
  }
  return defvalue;
}

// argumentStackKeyToDouble()
double argumentStackKeyToDouble(argStackType **argumentStackAddress, 
                                       char *match, double defvalue) {
  int arg;
  argStackType *argumentStack;
  argumentStack = *argumentStackAddress;
  for (arg = 0; arg < argumentStack->nargs - 1; arg++) {
    if ( argumentStack->marked[arg]=='F' && strcmp(argumentStack->arglist[arg], match)==0 ) {
      // Check if next arg already marked
      if (argumentStack->marked[arg+1]=='T') continue; 
      argumentStack->marked[arg] = 'T';
      arg++;
      argumentStack->marked[arg] = 'T';
      return atof(argumentStack->arglist[arg]);
    }
  }
  return defvalue;
}

// argumentStackKeyToFloat()
float argumentStackKeyToFloat(argStackType **argumentStackAddress,
                                       char *match, float defvalue) {
  int arg;
  argStackType *argumentStack;
  argumentStack = *argumentStackAddress;
  for (arg = 0; arg < argumentStack->nargs - 1; arg++) {
    if ( argumentStack->marked[arg]=='F' && strcmp(argumentStack->arglist[arg], match)==0 ) {
      // Check if next arg already marked
      if (argumentStack->marked[arg+1]=='T') continue;
      argumentStack->marked[arg] = 'T';
      arg++;
      argumentStack->marked[arg] = 'T';
      return atof(argumentStack->arglist[arg]);
    }
  }
  return defvalue;
}

// argumentStackKeyToString()
char* argumentStackKeyToString(argStackType **argumentStackAddress, 
                                      char *match, char* defvalue) {
  int arg;
  char *retvalue = NULL;
  argStackType *argumentStack;
  argumentStack = *argumentStackAddress;
  for (arg = 0; arg < argumentStack->nargs - 1; arg++) {
    if ( argumentStack->marked[arg]=='F' && strcmp(argumentStack->arglist[arg], match)==0 ) {
      // Check if next arg already marked
      if (argumentStack->marked[arg+1]=='T') continue; 
      argumentStack->marked[arg] = 'T';
      arg++;
      argumentStack->marked[arg] = 'T';
      retvalue = (char*) malloc( (strlen(argumentStack->arglist[arg])+1) * sizeof(char));
      strcpy(retvalue, argumentStack->arglist[arg]);
      return retvalue;
    }
  }
  if (defvalue!=NULL) {
    retvalue = (char*) malloc( (strlen(defvalue)+1) * sizeof(char));
    strcpy(retvalue, defvalue);
  }
  return retvalue;
}

// printArgumentStack()
void printArgumentStack(argStackType **argumentStackAddress) {
  int arg;
  argStackType *argumentStack;
  argumentStack = *argumentStackAddress;
  for (arg = 0; arg < argumentStack->nargs; arg++)
    fprintf(stdout,"Arg %i: %s [%c]\n",arg,argumentStack->arglist[arg],
            argumentStack->marked[arg]);
}

