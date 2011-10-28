#include "ArgList.h"
//ArgList.cpp
#include <cstring>
#include "CpptrajStdio.h"
#include <locale>
#include <iostream>
#include <sstream>
#include <stdexcept>
 
using namespace std;

// CONSTRUCTOR
ArgList::ArgList() {
  debug = 0;
}

// DESTRUCTOR
ArgList::~ArgList() {}

/* ArgList::SetDebug()
 * Set the arglist debug level.
 */
void ArgList::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0)
    mprintf("ArgList debug level set to %i\n",debug);
}

/* ArgList::SetList()
 * Separate input by the characters in separator and store as separate args.
 * This overwrites any existing args and completely resets the list.
 */
int ArgList::SetList(char *inputString, const char *separator) {
  string argument;

  if (inputString==NULL || separator==NULL) return 1;
  // Copy inputString to temp since it is destroyed by tokenize,
  // this allows const strings to be passed in.
  char *tempString = new char[strlen(inputString)+1];
  strcpy(tempString,inputString);
  
  // Free existing arglist
  arglist.clear();
  marked.clear();
  
  // Store inputString
  argline.assign(inputString);

  // Begin tokenization
  char *pch = strtok(tempString, separator);
  if (pch!=NULL) {
    while (pch!=NULL) {
      //if (debug>1) mprintf("getArgList:  Arg %i, Token [%s], ",nargs,pch);
      // If the argument is not quoted add it to the list
      if (pch[0]!='"') {
        argument.assign(pch);
        arglist.push_back(argument);

      // If the argument begins with a quote, place this and all subsequent
      // arguments ending with another quote into this argument
      } else {
        argument.assign(pch);
        // Check if this argument itself ends with a quote
        unsigned int argsize = argument.size();
        if (argsize == 1 || argument[argsize-1]!='"') {
          while (pch!=NULL) {
            argument.append(" ");
            pch=strtok(NULL," ");
            argument.append(pch);
            if (strchr(pch,'"')!=NULL) break;
          }
        }
        // Remove quotes from the argument
        for (string::iterator character = argument.begin();
                              character < argument.end();
                              character++)
          if (*character == '"') character = argument.erase(character);
        arglist.push_back(argument);
      }
      //if (debug>1) mprintf("Arglist[%i]= [%s]\n",nargs-1,arglist[nargs-1]);
      pch = strtok(NULL,separator);
    } // END while loop
    // Set up marked array
    marked.resize( arglist.size(), false );
  }
  // if (debug>0) mprintf("getArgList: Processed %i args\n",nargs);
  delete[] tempString;
  return 0;
}

/* ArgList::AddArg()
 * Add input to the argument list.
 */
void ArgList::AddArg(char *input) {
  string argument;
  // Dont store blank tokens
  if (input==NULL) return;
  if (input[0]=='\n') return;
  argument.assign(input);
  arglist.push_back(argument);
  argline.append(argument);
  argline.append(" ");
  marked.push_back(false);
}

/* ArgList::ResetMarked()
 * Reset all entries in the marked list to false. If noResetCmd is true
 * dont reset the first arg.
 * NOTE: Is this essential? If it is, is the bool essential?
 */
void ArgList::ResetMarked(bool noResetCmd) {
  unsigned int startarg = 0;
  if (noResetCmd) startarg = 1;
  for (unsigned int arg = startarg; arg < marked.size(); arg++)
    marked[arg]=false;
}

/* ArgList::MarkAll()
 * Set all entries in the marked list to true.
 */
void ArgList::MarkAll() {
  for (unsigned int arg = 0; arg < marked.size(); arg++)
    marked[arg]=true;
}

/* ArgList::CheckForMoreArgs()
 * Check if all arguments have been processed. If not print a warning along
 * with all unprocessed arguments.
 */
void ArgList::CheckForMoreArgs() {
  bool empty = true;
  string notmarked;
  
  for (unsigned int arg=0; arg < arglist.size(); arg++) {
    if (!marked[arg]) {
      empty=false;
      notmarked.append(arglist[arg] + " ");
    }
  }
  if (!empty)  
    mprintf("Warning: [%s] Not all arguments handled: [ %s]\n",arglist[0].c_str(),
            notmarked.c_str());
}

/* ArgList::PrintList()
 * Print out each arg on separate lines.
 */
void ArgList::PrintList() {
  unsigned int nargs = arglist.size();
  if (debug==0) {
    for (unsigned int arg = 0; arg < nargs; arg++) 
      mprintf("  %u: %s\n",arg,arglist[arg].c_str());
  } else {
    mprintf("ArgLine: %s\n",argline.c_str());
    for (unsigned int arg = 0; arg < nargs; arg++)
      mprintf("\tArg %u: %s (%i)\n",arg,arglist[arg].c_str(),(int)marked[arg]);
  }
}

/* ArgList::ArgLine()
 * Return the original argument string
 */
const char *ArgList::ArgLine() {
  return argline.c_str();
}        
        
/* ArgList::ArgAt()
 * Return arg at specified position.
 */
char *ArgList::ArgAt(int pos) {
  if (pos < 0 || pos >= (int) arglist.size()) return NULL;
  return (char*)arglist[pos].c_str();
}

/* ArgList::ArgIs()
 * Return true if arg at specified position matches input.
 */
bool ArgList::ArgIs(int pos, const char *input) {
  if (pos < 0 || pos >= (int) arglist.size()) return false;
  if (arglist[pos].compare( input )==0) return true;
  return false;
}

/* ArgList::Command()
 * Return the first argument and mark it, even if already marked.
 */
const char *ArgList::Command() {
  if (arglist.empty()) return NULL;
  marked[0]=true;
  return arglist[0].c_str();
}

/* ArgList::CommandIs()
 * Check if key is command, return true if so. Mark command no matter what.
 */
bool ArgList::CommandIs(const char *key) {
  if (arglist.empty()) return false;
  marked[0]=true;
  if (arglist[0].compare( key )==0) return true;
  return false;
}

/* ArgList::CommandIs()
 * Check the first nchar characters of key against command, return true if
 * they match. Mark command no matter what.
 */
/*bool ArgList::CommandIs(const char *key, size_t nchar) {
  if (arglist.empty()) return false;
  marked[0]=true;
  if (arglist[0].compare( 0, nchar, key )==0) return true;
  return false;
}*/

/* ArgList::getNextString()
 * Return the next unmarked string.
 * NOTE: The case back to char* is potentially dangerous if calling routine
 *       modifies the string in any way.
 */
char *ArgList::getNextString() {
  for (unsigned int arg = 0; arg < arglist.size(); arg++)
    if (!marked[arg]) {
      marked[arg]=true;
      return (char*)arglist[arg].c_str();
    }
  return NULL;
}

/* ArgList::getNextMask()
 * Return next unmarked Mask. A mask MUST include one of the following: 
 *   ':' residue
 *   '@' atom
 *   '*' everything
 *   NOTE: Disabling the following for now:
 *   '/' element
 *   '%' type
 */
char *ArgList::getNextMask() {
  for (unsigned int arg=0; arg < arglist.size(); arg++) {
    if (!marked[arg]) {
      const char *argmnt = arglist[arg].c_str();
      if ( strchr( argmnt, ':')!=NULL ||
           strchr( argmnt, '@')!=NULL ||
           strchr( argmnt, '*')!=NULL //||
           //strchr( argmnt, '/')!=NULL ||
           //strchr( argmnt, '%')!=NULL    
         )
      {
        marked[arg]=true;
        return (char*)argmnt;
      }
    }
  }
  return NULL;
}

/* Class: BadConversion
 * Runtime exception class for catching bad conversions from the 
 * convertToX routines.
 */
class BadConversion : public std::runtime_error {
public:
  BadConversion(std::string const &s)
    : std::runtime_error(s)
    { }
};

/* convertToInteger()
 * Convert the input string to an integer.
 */
inline int convertToInteger(string const &s) {
  istringstream iss(s);
  int i;
  if (!(iss >> i))
    throw BadConversion("convertToInteger(\"" + s + "\")");
  return i;
}

/* convertToDouble()
 * Convert the input string to a double.
 */
inline double convertToDouble(string const &s) {
  istringstream iss(s);
  double d;
  if (!(iss >> d))
    throw BadConversion("convertToDouble(\"" + s + "\")");
  return d;
}

/* validInteger()
 * Brief check that the passed in string begins with a digit
 * or '-'
 */
inline bool validInteger(string const &argument) {
  locale loc;
  if (isdigit(argument[0],loc) || argument[0]=='-') return true;
  return false;
}

/* validDouble()
 * Brief check that the passed in string begins with a digit,
 * '-', or '.'
 */
inline bool validDouble(string const &argument) {
  locale loc;
  if (isdigit(argument[0],loc) || argument[0]=='-' || argument[0]=='.' ) return true;
  return false;
}

/* ArgList::getNextInteger()
 * Convert next unmarked string to int and return, otherwise return def
 */
int ArgList::getNextInteger(int def) {
  locale loc;
  for (unsigned int arg=0; arg < arglist.size(); arg++)
    if (!marked[arg]) {
      // Check that first char is indeed an integer or '-', if not then continue
      if (validInteger(arglist[arg])) {
        int ival = convertToInteger(arglist[arg]);
        marked[arg]=true;
        return ival;
      }
    }
  return def;
}

/* ArgList::getNextDouble()
 * Convert next unmarked string to double and return, otherwise return def
 */
double ArgList::getNextDouble(double def) {
  locale loc;
  for (unsigned int arg=0; arg < arglist.size(); arg++)
    if (!marked[arg]) {
      // Check that first char is indeed a digit, '.', or '-', if not then continue
      if (validDouble(arglist[arg])) {
        double dval = convertToDouble(arglist[arg]);
        marked[arg]=true;
        return dval;
      }
    }
  return def;
}

/* ArgList::getKeyString()
 * Search for unmarked key in arglist, return if found, otherwise return def
 */
char *ArgList::getKeyString(const char *key, char *def) {
  unsigned int nargs = arglist.size() - 1;
  for (unsigned int arg=0; arg < nargs; arg++)
    if (!marked[arg]) {
      if (arglist[arg].compare(key)==0) { 
        marked[arg]=true;
        arg++;
        marked[arg]=true;
        return (char*)arglist[arg].c_str();
      }
    }
  return def;
}

/* ArgList::getKeyIndex()
 * Search for unmarked key in arglist, return its position in the ArgList
 * if found, otherwise return -1
 */
/*
int ArgList::getKeyIndex(char *key) {
  for (unsigned int arg=0; arg < arglist.size(); arg++) {
    if (arglist[arg].compare(key)==0) return arg;
  }
  return -1;
}
*/

/* ArgList::getKeyInt()
 * Search for unmarked key in arglist, return if found, otherwise returh def
 */
int ArgList::getKeyInt(const char *key, int def) {
  unsigned int nargs = arglist.size() - 1;
  for (unsigned int arg=0; arg < nargs; arg++)
    if (!marked[arg]) {
      if (arglist[arg].compare(key)==0) {
        if (validInteger(arglist[arg+1])) {
          marked[arg]=true;
          arg++;
          int ival = convertToInteger(arglist[arg]);
          marked[arg]=true;
          return ival;
        }
      }
    }
  return def;
}

/* ArgList::getKeyDouble()
 * Search for unmarked key in arglist, return if found, otherwise return def
 */
double ArgList::getKeyDouble(const char *key, double def) {
  unsigned int nargs = arglist.size() - 1;
  for (unsigned int arg=0; arg < nargs; arg++)
    if (!marked[arg]) {
      if (arglist[arg].compare(key)==0) {
        if (validDouble(arglist[arg+1])) {
          marked[arg]=true;
          arg++;
          double dval = convertToDouble(arglist[arg]);
          marked[arg]=true;
          return dval;
        }
      }
    }
  return def;
}

/* ArgList::hasKey()
 * Return true if key is found, false if not. Mark key if found.
 */
bool ArgList::hasKey(const char *key) {
  for (unsigned int arg = 0; arg < arglist.size(); arg++) 
    if (!marked[arg]) {
      if (arglist[arg].compare(key)==0) {
        marked[arg]=true;
        return true;
      }
    }
  return false;
}

/* ArgList::Contains()
 * Like hasKey(), but key is not marked. Return true if key is found, 
 * false if not.
 * NOTE: Should this be ignoring previously marked strings?
 */
bool ArgList::Contains(const char *key) {
  for (unsigned int arg = 0; arg < arglist.size(); arg++) 
    if (!marked[arg]) {
      if (arglist[arg].compare(key)==0) {
        return true;
      }
    }
  return false;
}


