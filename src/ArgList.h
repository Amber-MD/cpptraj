#ifndef INC_ARGLIST_H
#define INC_ARGLIST_H
/// Class: ArgList
/// Test of ArgList using vector and string STL classes.
/// Hold a list of string arguments, can be set from a delimited list
/// with SetList or arguments can be added individually with AddArg. 
/// Arguments can be accessed with the various getX routines,
/// where X is specific to certain types, e.g. getNextDouble returns
/// the next double, getNextMask returns an atom mask expression (i.e.
/// it has :, @, % characters etc).
#include <vector>
#include <string>
class ArgList {
    std::vector<std::string> arglist;
    std::string argline;
    std::vector<bool> marked;
    int debug;
  public:
    ArgList();
    ~ArgList();
    void SetDebug(int);

    int SetList(char *, const char *);
    void AddArg(char*);
    void ResetMarked(bool);
    void MarkAll();
    void CheckForMoreArgs();
    
    void PrintList();

    const char *ArgLine();
    char *ArgAt(int);
    bool ArgIs(int,const char*);
    const char *Command();
    bool CommandIs(const char*);
    //bool CommandIs(const char*,size_t);
    int Nargs() { return (int)arglist.size(); }

    char *getNextString();
    char *getNextMask();
    int getNextInteger(int);
    double getNextDouble(double);
    char *getKeyString(const char *, char*);
    //int getKeyIndex(char*); 
    int getKeyInt(const char *, int);
    double getKeyDouble(const char*, double);
    bool hasKey(const char*);
    bool Contains(const char*);
};
#endif
