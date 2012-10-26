#include "ArgList.h"
//ArgList.cpp
#include <cstring>
#include <locale>
#include <stdexcept> // ArgList[]
#include "CpptrajStdio.h"
#include "StringRoutines.h"
 
using namespace std;

// CONSTRUCTOR
ArgList::ArgList() :
  debug(0)
{}

// CONSTRUCTOR - Take string and convert to args delimited by space
ArgList::ArgList(const char* argIn) :
  debug(0)
{
  SetList(argIn, " ");
}

// CONSTRUCTOR - Take string and convert to args delimted by separators
ArgList::ArgList(std::string const& expression, const char *separators) :
  debug(0)
{
  if (!expression.empty())
    SetList(expression.c_str(), separators);
}

// COPY CONSTRUCTOR
ArgList::ArgList(const ArgList &rhs) {
  arglist = rhs.arglist;
  argline = rhs.argline;
  marked = rhs.marked;
  debug = rhs.debug;
}

// ArgList::operator=()
ArgList &ArgList::operator=(const ArgList &rhs) {
  if (&rhs==this) return *this;
  // Allocate and copy. Assignment ops should automatically deallocate.
  arglist = rhs.arglist;
  argline = rhs.argline;
  marked = rhs.marked;
  debug = rhs.debug;
  return *this;
}

// ArgList::operator[]
std::string const& ArgList::operator[](int idx) {
  if (idx < 0 || idx >= (int)arglist.size())
    throw std::out_of_range("ArgList[]");
  return arglist[idx];
}

// ArgList::SetDebug()
void ArgList::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0)
    mprintf("ArgList debug level set to %i\n",debug);
}

// ArgList::SetList()
int ArgList::SetList(const char *inputString) {
  return SetList(inputString, " ");
}

// ArgList::SetList()
int ArgList::SetList(std::string const& inputString, const char* separator) {
  if (inputString.empty()) return 1;
  return (SetList(inputString.c_str(), separator));
}

// ArgList::SetList()
/** Separate input by the characters in separator and store as separate args.
  * This overwrites any existing args and completely resets the list.
  * \param inputString NULL terminated string to be converted to arguments
  * \param separator string containing characters used to separate arguments
  * \return 0 if arglist successfully set up, 1 if not.
  */
int ArgList::SetList(const char *inputString, const char *separator) {
  string argument;
  char quotechar;

  if (inputString==NULL || separator==NULL) return 1;
  // Copy inputString to temp since it is destroyed by tokenize,
  // this allows const strings to be passed in.
  size_t inputStringSize = strlen(inputString);
  char *tempString = new char[ inputStringSize+1 ];
  if (inputStringSize < 1) return 1;
  strcpy(tempString,inputString);
  // Remove newline char from tempString if present
  if ( tempString[ inputStringSize - 1 ] == '\n' )
    tempString[ inputStringSize - 1 ] = '\0';
  
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
        unsigned int argsize = argument.size();
        // Check for blank quote token ("")
        if (argsize != 2 || argument[1] != quotechar) {
          // Check if this argument itself ends with a quote
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
void ArgList::AddArg(const char *input) {
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
      mprintf("  %u: %s\n",arg+1,arglist[arg].c_str());
  } else {
    mprintf("ArgLine: %s\n",argline.c_str());
    for (unsigned int arg = 0; arg < nargs; arg++)
      mprintf("\tArg %u: %s (%i)\n",arg+1,arglist[arg].c_str(),(int)marked[arg]);
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
const char *ArgList::ArgAt(int pos) {
  if (pos < 0 || pos >= (int) arglist.size()) return NULL;
  return arglist[pos].c_str();
}

void ArgList::RemoveFirstArg() {
  if (arglist.empty()) return;
  arglist.erase( arglist.begin() );
  marked.erase( marked.begin() );
}

// ArgList::Command()
/* \return pointer to the first argument
 */
const char *ArgList::Command() const {
  if (arglist.empty()) return NULL;
  return arglist[0].c_str();
}

// ArgList::CommandIs()
/** \param key Key to check first argument against
  * \return true if first argument matches key
  */
bool ArgList::CommandIs(const char *key) const {
  if (arglist.empty()) return false;
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
ArgList::ConstArg ArgList::getNextString() {
  for (unsigned int arg = 0; arg < arglist.size(); ++arg)
    if (!marked[arg]) {
      marked[arg]=true;
      return arglist[arg].c_str();
    }
  return NULL;
}

// ArgList::GetStringNext()
std::string ArgList::GetStringNext() {
  std::string emptystring;
  for (unsigned int arg = 0; arg < arglist.size(); ++arg)
    if (!marked[arg]) {
      marked[arg]=true;
      return arglist[arg];
    }
  return emptystring;
}

// ArgList::getNextMask()
/** Return next unmarked Mask. A mask MUST include one of the following: 
  *   ':' residue
  *   '@' atom
  *   '/' element
  *   '%' type
  *   '*' everything
  * \return the next unmarked atom mask expression
  */
ArgList::ConstArg ArgList::getNextMask() {
  for (unsigned int arg=0; arg < arglist.size(); ++arg) {
    if (!marked[arg]) {
      size_t found = arglist[arg].find_first_of(":@*/%");
      if (found != std::string::npos) {
        marked[arg]=true;
        return arglist[arg].c_str();
      }
    }
  }
  return NULL;
}

// ArgList::GetMaskNext()
std::string ArgList::GetMaskNext() {
  for (unsigned int arg = 0; arg < arglist.size(); ++arg) {
    if (!marked[arg]) {
      size_t found = arglist[arg].find_first_of(":@*/%");
      if (found != std::string::npos) {
        marked[arg] = true;
        return arglist[arg];
      }
    }
  }
  return std::string();
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
  * as a string if found, otherwise return 0.
  * \param key string to search for
  */
ArgList::ConstArg ArgList::getKeyString(const char *key) {
  unsigned int nargs = arglist.size() - 1;
  for (unsigned int arg=0; arg < nargs; ++arg)
    if (!marked[arg]) {
      if (arglist[arg].compare(key)==0) { 
        marked[arg++]=true;
        marked[arg]=true;
        return arglist[arg].c_str();
      }
    }
  return NULL;
}

// ArgList::GetStringKey()
std::string ArgList::GetStringKey(const char *key) {
  std::string empty;
  unsigned int nargs = arglist.size() - 1;
  for (unsigned int arg=0; arg < nargs; arg++)
    if (!marked[arg]) {
      if (arglist[arg].compare(key)==0) {
        marked[arg]=true;
        arg++;
        marked[arg]=true;
        return arglist[arg];
      }
    }
  return empty;
}

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

// ArgList::ArgToDouble()
double ArgList::ArgToDouble(int pos) {
  if (pos < 0 || pos >= (int)arglist.size()) {
    mprinterr("Internal Error: ArgList::ArgToDouble: position out of range.\n");
    return 0;
  }
  if (!validDouble(arglist[pos])) {
    mprinterr("Error: Arg %s is not a valid double.\n",arglist[pos].c_str());
    return 0;
  }
  return convertToDouble(arglist[pos]);
}

// ArgList::ArgToInteger()
int ArgList::ArgToInteger(int pos) {
  if (pos < 0 || pos >= (int)arglist.size()) {
    mprinterr("Internal Error: ArgList::ArgToInteger: position out of range.\n");
    return 0;
  }
  if (!validInteger(arglist[pos])) {
    mprinterr("Error: Arg %s is not a valid integer.\n",arglist[pos].c_str());
    return 0;
  }
  return convertToInteger(arglist[pos]);
}

/** Return true as soon as an unmarked arg is encountered. */
bool ArgList::ArgsRemain() {
  for (std::vector<bool>::iterator mark = marked.begin();
                                   mark != marked.end(); ++mark)
  {  
    if ( !(*mark) ) return true;
  }
  return false;
}
