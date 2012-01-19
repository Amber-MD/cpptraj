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

// ArgList::SetDebug()
void ArgList::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0)
    mprintf("ArgList debug level set to %i\n",debug);
}

// ArgList::SetList()
/** Separate input by the characters in separator and store as separate args.
  * This overwrites any existing args and completely resets the list.
  * \param inputString NULL terminated string to be converted to arguments
  * \param separator string containing characters used to separate arguments
  * \return 0 if arglist successfully set up, 1 if not.
  */
int ArgList::SetList(char *inputString, const char *separator) {
  string argument;
  char quotechar;

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
      if      (pch[0]=='"' ) quotechar='"';
      else if (pch[0]=='\'') quotechar='\'';
      else quotechar=' ';
      if (quotechar==' ') {
        argument.assign(pch);
        arglist.push_back(argument);

      // If the argument begins with a quote, place this and all subsequent
      // arguments ending with another quote into this argument
      } else {
        argument.assign(pch);
        // Check if this argument itself ends with a quote
        unsigned int argsize = argument.size();
        if (argsize == 1 || argument[argsize-1]!=quotechar) {
          while (pch!=NULL) {
            argument.append(" ");
            pch=strtok(NULL," ");
            // If pch is NULL at this point there was no closing quote.
            if (pch==NULL) {
              mprintf("\tWarning: argument missing closing quote [%c]\n",quotechar);
              break;
            }
            argument.append(pch);
            if (strchr(pch,quotechar)!=NULL) break;
          }
        }
        // Remove quotes from the argument
        for (string::iterator character = argument.begin();
                              character < argument.end();
                              character++)
          if (*character == quotechar) character = argument.erase(character);
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

// ArgList::AddArg()
/** \param input NULL terminated string to add to argument list
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

// ArgList::ResetMarked()
/** \param noResetCmd If true dont reset the marked status of the first arg.
  */
// NOTE: Is this essential? If it is, is the bool essential?
void ArgList::ResetMarked(bool noResetCmd) {
  unsigned int startarg = 0;
  if (noResetCmd) startarg = 1;
  for (unsigned int arg = startarg; arg < marked.size(); arg++)
    marked[arg]=false;
}

// ArgList::MarkAll()
void ArgList::MarkAll() {
  for (unsigned int arg = 0; arg < marked.size(); arg++)
    marked[arg]=true;
}

// ArgList::MarkArg()
void ArgList::MarkArg(int arg) {
  if (arg < 0 || arg >= (int) marked.size()) return;
  marked[arg]=true;
}

// ArgList::CheckForMoreArgs()
/** Check if all arguments have been processed. If not print a warning along
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

// ArgList::PrintList()
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

// ArgList::ArgLine()
const char *ArgList::ArgLine() {
  return argline.c_str();
}        
        
// ArgList::ArgAt()
/** \param pos argument position
  * \return pointer to the argument at pos or NULL if pos is out of bounds.
  */
char *ArgList::ArgAt(int pos) {
  if (pos < 0 || pos >= (int) arglist.size()) return NULL;
  return (char*)arglist[pos].c_str();
}

// ArgList::ArgIs()
/** \param pos argument position
  * \param input Key to check arguments against.
  * \return true if argument at pos is input
  */
bool ArgList::ArgIs(int pos, const char *input) {
  if (pos < 0 || pos >= (int) arglist.size()) return false;
  if (arglist[pos].compare( input )==0) return true;
  return false;
}

// ArgList::Command()
/** Return the first argument and mark it, even if already marked.
 * \return pointer to the first argument
 */
const char *ArgList::Command() {
  if (arglist.empty()) return NULL;
  marked[0]=true;
  return arglist[0].c_str();
}

// ArgList::CommandIs()
/** Check if key matches the first argument. Mark command no matter what.
  * \param key Key to check first argument against
  * \return true if first argument matches key
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

// ArgList::getNextString()
/** \return the next unmarked string.
  */
// NOTE: The case back to char* is potentially dangerous if calling routine
//        modifies the string in any way.
char *ArgList::getNextString() {
  for (unsigned int arg = 0; arg < arglist.size(); arg++)
    if (!marked[arg]) {
      marked[arg]=true;
      return (char*)arglist[arg].c_str();
    }
  return NULL;
}

// ArgList::getNextMask()
/** Return next unmarked Mask. A mask MUST include one of the following: 
  *   ':' residue
  *   '@' atom
  *   '*' everything
  * \return the next unmarked atom mask expression
  */
// NOTE: Disabling the following for now:
// '/' element
// '%' type
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

// ArgList::getNextTag()
/** Return the next unmarked tag. A tag is defined as a character string
  * bounded by brackets, e.g. [tag].
  */
string ArgList::getNextTag() {
  string emptystring;
  for (unsigned int arg = 0; arg < arglist.size(); arg++) {
    if (!marked[arg]) {
      string::reverse_iterator lastchar  = arglist[arg].rbegin();
      string::iterator         firstchar = arglist[arg].begin();
      if (*firstchar=='[' && *lastchar==']') {
        marked[arg]==true;
        return arglist[arg];
      }
    }
  }
  return emptystring;
}

/*! \class: BadConversion
    \brief Runtime exception for catching bad conversions from the convertToX routines.
  */
class BadConversion : public std::runtime_error {
public:
  BadConversion(std::string const &s)
    : std::runtime_error(s)
    { }
};

// convertToInteger()
/// Convert the input string to an integer.
inline int convertToInteger(string const &s) {
  istringstream iss(s);
  long int i;
  if (!(iss >> i))
    throw BadConversion("convertToInteger(\"" + s + "\")");
  return (int)i;
}

// convertToDouble()
/// Convert the input string to a double.
inline double convertToDouble(string const &s) {
  istringstream iss(s);
  double d;
  if (!(iss >> d))
    throw BadConversion("convertToDouble(\"" + s + "\")");
  return d;
}

// validInteger()
/// Brief check that the passed in string begins with a digit or '-'
inline bool validInteger(string const &argument) {
  locale loc;
  if (isdigit(argument[0],loc) || argument[0]=='-') return true;
  return false;
}

// validDouble()
/// Brief check that the passed in string begins with a digit, '-', or '.'
inline bool validDouble(string const &argument) {
  locale loc;
  if (isdigit(argument[0],loc) || argument[0]=='-' || argument[0]=='.' ) return true;
  return false;
}

// ArgList::getNextInteger()
/** \param def Value to return if no integer args found
  * \return Next unmarked integer argument or def
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

// ArgList::getNextDouble()
/** \param def Value to return if no double args found
  * \return Next unmarked double argument or def
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

// ArgList::getKeyString()
/** Search the argument list for key, return the argument following key
  * as a string if found, otherwise return def.
  * \param key string to search for
  * \param def Value to return if key not found.
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

// ArgList::getKeyInt()
/** Search the argument list for key, return the argument following key
  * as an integer if found, otherwise return def.
  * \param key string to search for
  * \param def Value to return if key not found.
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

// ArgList::getKeyDouble()
/** Search the argument list for key, return the argument following key
  * as a double if found, otherwise return def.
  * \param key string to search for
  * \param def Value to return if key not found.
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

// ArgList::hasKey()
/** Search the argument list for key, mark and return true if found.
  * \param key string to search for
  * \return true if key is found, false if not.
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

// ArgList::Contains()
/** \param key string to search for
  * \return true if key is found, false if not.
  */
// NOTE: Should this be ignoring previously marked strings?
bool ArgList::Contains(const char *key) {
  for (unsigned int arg = 0; arg < arglist.size(); arg++) 
    if (!marked[arg]) {
      if (arglist[arg].compare(key)==0) {
        return true;
      }
    }
  return false;
}


