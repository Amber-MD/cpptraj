#ifndef INC_ARGLIST_H
#define INC_ARGLIST_H
#include <vector>
#include <string>
// Class: ArgList
/// Hold a list of string arguments and keeps track of their usage.
/** Can be set from an input line using SetList(), with arguments separated 
  * by a specified delimiter, or arguments can be added one-by-one with AddArg.
  * Arguments can be accessed with the various getX routines,
  * where X is specific to certain types, e.g. getNextDouble returns
  * the next double, getNextMask returns an atom mask expression (i.e.
  * it has :, @, % characters etc). All of the getX routines (along with
  * the hasKey routine) mark the argument they access as used, so that
  * subsequent calls with these functions will not return the same
  * argument over and over. For string arguments two versions currently
  * exist, lowercase 'get' and uppercase 'Get'. The 'get' routines
  * return a const char* (type defined as ConstArg), while the 'Get' routines
  * return a string object. The former routines exist both as legacy code
  * and potentially for performance reasons (since no copy of a string
  * object has to be made), although I haven't benchmarked to see just
  * how "fast" they are). It should only be used in local scope, i.e.
  * to hold a variable in Action::init() only; anything outside that
  * should use the 'Get' routine.
  */
class ArgList {
  public:
    /// Return type definition for 'getX' string routines.
    typedef const char* ConstArg;
    // Construction/Assignment/Access
    ArgList();
    ArgList(const char*);
    ArgList(std::string const&, const char*);
    ArgList(const ArgList&);
    ArgList& operator=(const ArgList &);
    std::string const& operator[](int);
    // Iterators
    typedef std::vector<std::string>::const_iterator const_iterator;
    const_iterator begin() { return arglist.begin(); }
    const_iterator end()   { return arglist.end();   }
    /// Set the debug level
    void SetDebug(int);
    /// Set up argument list from string using space as a separator.
    int SetList(const char *);
    /// Set up argument list from string and given separators
    int SetList(const char *, const char *); // TODO: Obsolete
    /// Set up argument list from string and given separators
    int SetList(std::string const&, const char *);
    /// Add argument to the list
    void AddArg(const char*);
    /// Mark given argument
    void MarkArg(int);
    /// Print a warning if not all arguments are marked
    void CheckForMoreArgs();
    /// Print the argument list
    void PrintList();
    /// Return the argument string
    const char *ArgLine();
    /// Return the argument at given position
    const char* ArgAt(int);
    void RemoveFirstArg();
    /// Return the first argument
    const char *Command() const;
    /// Return true if the first argument matches key
    bool CommandIs(const char*) const;
    //bool CommandIs(const char*,size_t);
    /// Return the number of arguments
    int Nargs()  { return (int)arglist.size(); }
    bool empty() { return arglist.empty();     }
    /// Return const char* to next unmarked string.
    ConstArg getNextString();
    /// Return the next unmarked string
    std::string GetStringNext();
    /// Return the next unmarked mask
    ConstArg getNextMask();
    /// Return the next unmarked mask
    std::string GetMaskNext();
    /// Return the next unmarked tag
    std::string getNextTag();
    /// Return the next unmarked integer
    int getNextInteger(int);
    /// Return the next unmarked double
    double getNextDouble(double);
    /// Return const char* to string following the given key.
    ConstArg getKeyString(const char *);
    /// Return the string following the given key
    std::string GetStringKey(const char *);
    /// Return the integer following the given key 
    int getKeyInt(const char *, int);
    /// Return the double following the given key
    double getKeyDouble(const char*, double);
    /// Return true if the key is present in the list
    bool hasKey(const char*);
    /// Return true if they key is in the list but do not mark.
    bool Contains(const char*);
    /// Convert arg at position to double.
    double ArgToDouble(int);
    /// Convert arg at position to integer.
    int ArgToInteger(int);
    /// True if any unmarked args remain.
    bool ArgsRemain();
  private:
    /// List of arguments
    std::vector<std::string> arglist;
    /// The original argument string (complete list)
    std::string argline;
    /// Mark which arguments have been used
    std::vector<bool> marked;
    /// debug level
    int debug;
};
#endif
