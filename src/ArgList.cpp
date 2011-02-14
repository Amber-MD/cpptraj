// ArgList
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <string> // Used in constructor to detect quoted args
#include "ArgList.h"
#include "CpptrajStdio.h"
// NOTE: Place checks on memory

// CONSTRUCTOR: Empty Arg List
ArgList::ArgList() {
  nargs=0;
  arglist=NULL;
  marked=NULL;
  argline=NULL;
}

// CONSTRUCTOR
// Separate input by the characters in separator and store as separate args
ArgList::ArgList(char *input, const char *separator) {
  char *pch;
  std::string quotedArg;
  int debug;
  size_t inputSize;

  arglist=NULL; marked=NULL; argline=NULL; nargs=0; debug=0;

  inputSize = strlen(input);
  // Replace any trailing newline char from input with NULL
  //nargs=strlen(input);
  if (inputSize>0) {
    if (input[inputSize-1]=='\n') input[inputSize-1]='\0';
  }
  if (debug>3) 
    mprintf("getArgList: Setting up arg list for [%s] with separator [%s]\n",
            input, separator);

  // Store original argument line
  argline=(char*) malloc( (inputSize+1) * sizeof(char) );
  strcpy(argline, input);

  pch=strtok(input,separator);
  if (pch!=NULL) {

    while (pch!=NULL) {
      if (debug>3) mprintf("getArgList:  Arg %i, Token [%s], ",nargs,pch);
      if ( pch[0]!='"' ) 
        Add(pch);
      else {
        // If the argument begins with a quote, place this and all subsequent
        // arguments ending with another quote into the same argument.
        quotedArg.clear();
        quotedArg.assign(pch);
        // Check if this argument itself ends with a quote
        if (quotedArg.size() == 1 || quotedArg[quotedArg.size()-1]!='"') {
          while (pch!=NULL) {
            quotedArg.append(" ");
            pch=strtok(NULL," ");
            quotedArg.append(pch);
            if (strchr(pch,'"')!=NULL) break;
          }
        }
        // Remove the quotes from the argument
        for (std::string::iterator it=quotedArg.begin(); it < quotedArg.end(); it++)
          if (*it=='"') quotedArg.erase(it);  
        Add((char*)quotedArg.c_str());
      }
      if (debug>3) mprintf("Arglist[%i]= [%s]\n",nargs-1,arglist[nargs-1]);
      pch=strtok(NULL,separator);
    }
    // Setup marked array
    Reset();
  }
  if (debug>3) mprintf("getArgList: Processed %i args\n",nargs);
}

// DESTRUCTOR
ArgList::~ArgList() {
  int i;

  if (arglist!=NULL){

    for (i=0; i<nargs; i++)
      if (arglist[i]!=NULL) free(arglist[i]);
    free(arglist);
  }
  if (marked!=NULL) free(marked);
  if (argline!=NULL) free(argline);
}


/* 
 * ArgList::Add()
 * Add input to the argument list.
 */
void ArgList::Add(char *input) {
  // Dont store blank tokens
  if (input==NULL) return;
  if (input[0]!='\n') {
    arglist=(char**) realloc(arglist,(nargs+1)*sizeof(char*));
    arglist[nargs]=(char*) malloc( (strlen(input)+1) * sizeof(char) );
    strcpy(arglist[nargs],input);
    nargs++;
  }
}

/*
 * ArgList::Copy()
 * Return a copy of this arglist
 */
ArgList *ArgList::Copy() {
  ArgList *temp; 
  int i;

  temp=new ArgList();

  for (i=0; i<nargs; i++) 
    temp->Add(arglist[i]);

  if (argline!=NULL) {
    temp->argline=(char*) malloc( (strlen(argline)+1) * sizeof(char));
    strcpy(temp->argline,argline);
  }

  temp->Reset();

  return temp;
}

/*
 * ArgList::print()
 * Print out each arg on separate lines
 */
void ArgList::print() {
  int i;
  for (i=0; i<nargs; i++) 
    mprintf("  %i: %s\n",i,arglist[i]);
    //mprintf("  ArgList[%i]=%s\n",i,arglist[i]);
}

/*
 * ArgList::ArgLine()
 * Return the original argument string
 */
char *ArgList::ArgLine() {
  return argline;
}

/*
 * ArgList::Command()
 * Check the first arg for command
 * Mark and return. Return even if marked.
 */
char *ArgList::Command() {

  if (nargs==0) return NULL;
  marked[0]='T';
  return arglist[0];
}

/* ArgList::CommandIs()
 * Check if key is command, return 1 if true
 */
int ArgList::CommandIs(const char *key) {
  if (strcmp(this->Command(),key)==0) return 1;
  return 0;
}

/*
 * ArgList::CommandIs()
 * Check if first N chars of Command are key, return 1 if true
 */
int ArgList::CommandIs(const char *key, int N) {
  if (strncmp(this->Command(),key,N)==0) return 1;
  return 0;
}

/*
 * ArgList::getNextString()
 * Return next unmarked string
 */
char *ArgList::getNextString() {
  int i;
  for (i=0; i<nargs; i++)
    if (marked[i]!='T') {
      marked[i]='T';
      return arglist[i];
    }
  return NULL;
}

/*
 * ArgList::CheckForMoreArgs()
 * Check if all arguments have been processed. If not print a warning along
 * with all unprocessed arguments.
 */
void ArgList::CheckForMoreArgs() {
  int i;
  bool empty;

  empty=true;
  for (i=0; i<nargs; i++) {
    if (marked[i]=='F') {
      empty=false;
      break;
    }
  }
  if (!empty) {
    mprintf("Warning: [%s] Not all arguments handled: [ ",arglist[0]);
    for (i=0; i<nargs; i++) {
      if (marked[i]=='F') mprintf("%s ",arglist[i]);
    }
    mprintf(" ]\n");
  }
}

/* 
 * ArgList::getNextMask()
 * Return next unmarked Mask. A mask MUST include one of the following: 
 *   ':' residue
 *   '@' atom
 *   '*' everything
 *   NOTE: Disabling the following for now:
 *   '/' element
 *   '%' type
 */
char *ArgList::getNextMask() {
  int i;
  for (i=0; i<nargs; i++) {
    if (marked[i]!='T') {
      if ( strchr( arglist[i], ':')!=NULL ||
           strchr( arglist[i], '@')!=NULL ||
           strchr( arglist[i], '*')!=NULL //||
           //strchr( arglist[i], '/')!=NULL ||
           //strchr( arglist[i], '%')!=NULL    
         ) 
      {
        marked[i]='T';
        return arglist[i];
      }
    }
  }
  return NULL;
}

/*
 * ArgList::getNextInteger()
 * Convert next unmarked string to int and return, otherwise return def
 */
int ArgList::getNextInteger(int def) {
  int i;
  for (i=0; i<nargs; i++)
    if (marked[i]!='T') {
      // Check that first char is indeed an integer - if not continue
      if (!isdigit(arglist[i][0])) {
        //mprintf("WARNING: Getting integer from arg (%s) that is not digit!\n",arglist[i]);
        continue;
      }
      marked[i]='T';
      return atoi(arglist[i]);
    }
  return def;
}

/*
 * ArgList::getKeyString()
 * Search for unmarked key in arglist, return if found, otherwise return def
 */
char *ArgList::getKeyString(const char *key, char *def) {
  int i;

  for (i=0; i<nargs-1; i++)
    if (marked[i]!='T' && strcmp(key,arglist[i])==0) {
      marked[i]='T';
      i++;
      marked[i]='T';
      return arglist[i];
    }
  return def;
}

/*
 * ArgList::getKeyIndex()
 * Search for unmarked key in arglist, return its position in the ArgList
 * if found, otherwise return -1
 */
int ArgList::getKeyIndex(char *key) {
  int i;

  for (i=0; i<nargs; i++) {
    if (strcmp(key,arglist[i])==0) return i;
  }
  return -1;
}

/*
 * ArgList::getKeyInt()
 * Search for unmarked key in arglist, return if found, otherwise returh def
 */
int ArgList::getKeyInt(const char *key, int def) {
  int i;

  for (i=0; i<nargs-1; i++)
    if (marked[i]!='T' && strcmp(key,arglist[i])==0) {
      marked[i]='T'; 
      i++;
      marked[i]='T';
      // Brief check that first char is indeed an integer
      if (!isdigit(arglist[i][0])) {
        mprintf("WARNING: Getting integer from arg (%s) that is not digit!\n",arglist[i]);
      }
      return atoi(arglist[i]);
    }
  return def;
}

/*
 * ArgList::getKeyDouble()
 * Search for unmarked key in arglist, return if found, otherwise returh def
 */
double ArgList::getKeyDouble(const char *key, double def) {
  int i;

  for (i=0; i<nargs-1; i++)
    if (marked[i]!='T' && strcmp(key,arglist[i])==0) {
      marked[i]='T';
      i++;
      marked[i]='T';
      // Brief check that first char is indeed a digit or . 
      if (!isdigit(arglist[i][0]) && arglist[i][0]!='.') {
        mprintf("WARNING: getKeyDouble: arg (%s) does not appear to be a number!\n",
                arglist[i]);
      }
      return atof(arglist[i]);
    }
  return def;
}

/*
 * ArgList::hasKey()
 * Return 1 if key is found, 0 if not
 */
int ArgList::hasKey(const char *key) {
  int i;

  for (i=0; i<nargs; i++) {
    if (marked[i]!='T' && strcmp(key,arglist[i])==0) {
      marked[i]='T';
      return 1;
    }
  }
  return 0;
}

/*
 * ArgList::Reset()
 * Reset marked except for the Command (arg 0). Allocate memory for marked
 * if not already done.
 */
void ArgList::Reset() {
  if (nargs==0) return;
  if (marked==NULL) 
    marked=(char*) malloc(nargs*sizeof(char));
  if (marked!=NULL) {
    memset(marked,'F',nargs);
    marked[0]='T';
  }
}

/*
 * ArgList::ResetAll()
 * Reset all arguments including Command (arg 0)
 */
void ArgList::ResetAll() {
  if (nargs==0) return;
  if (marked==NULL)
    marked=(char*) malloc(nargs*sizeof(char));
  if (marked!=NULL) 
    memset(marked,'F',nargs);
}
    
