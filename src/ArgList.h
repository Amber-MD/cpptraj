#ifndef INC_ARGLIST_H
#define INC_ARGLIST_H
#include <vector>
#include <string>
// Class: ArgList
/// Hold a list of string arguments
/** Can be set from a delimited list
  * with SetList or arguments can be added individually with AddArg. 
  * Arguments can be accessed with the various getX routines,
  * where X is specific to certain types, e.g. getNextDouble returns
  * the next double, getNextMask returns an atom mask expression (i.e.
  * it has :, @, % characters etc). All of the getX routines (along with
  * the hasKey routine) mark the argument they access as used, so that
  * subsequent calls with these functions will not return the same
  * argument over and over.
  */
class ArgList {
    /// List of arguments
    std::vector<std::string> arglist;
    /// The original argument string (complete list)
    std::string argline;
    /// Mark which arguments have been used
    std::vector<bool> marked;
    /// debug level
    int debug;
  public:
    ArgList();
    ~ArgList();
    ArgList(const ArgList &);
    ArgList&operator=(const ArgList &);
    std::string& operator[](int);
    /// Set the debug level
    void SetDebug(int);
    /// Set up argument list from string and given separators
    int SetList(char *, const char *);
    /// Add argument to the list
    void AddArg(char*);
    /// Unmark all arguments
    void ResetMarked(bool);
    /// Mark all arguments
    void MarkAll();
    /// Mark given argument
    void MarkArg(int);
    /// Print a warning if not all arguments are marked
    void CheckForMoreArgs();
    /// Print the argument list
    void PrintList();
    /// Return the argument string
    const char *ArgLine();
    /// Return the argument at given position
    char *ArgAt(int);
    /// Return true if the argument at given position matches key
    bool ArgIs(int,const char*);
    /// Return the first argument
    const char *Command();
    /// Return true if the first argument matches key
    bool CommandIs(const char*);
    //bool CommandIs(const char*,size_t);
    /// Return the number of arguments
    int Nargs() { return (int)arglist.size(); }
    /// Return the next unmarked string
    char *getNextString();
    /// Return the next unmarked mask
    char *getNextMask();
    /// Return the next unmarked tag
    std::string getNextTag();
    /// Return the next unmarked integer
    int getNextInteger(int);
    /// Return the next unmarked double
    double getNextDouble(double);
    /// Return the string following the given key
    char *getKeyString(const char *, char*);
    //int getKeyIndex(char*);
    /// Return the integer following the given key 
    int getKeyInt(const char *, int);
    /// Return the double following the given key
    double getKeyDouble(const char*, double);
    /// Return the comma-separated arg following the given key
    ArgList getKeyArgList(const char *);
    /// Return true if the key is present in the list
    bool hasKey(const char*);
    /// Return true if they key is in the list but do not mark.
    bool Contains(const char*);
    /// Convert arg at position to double.
    double ArgToDouble(int);
    /// Convert arg at position to integer.
    int ArgToInteger(int);
};
#endif
