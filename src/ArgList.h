#ifndef INC_ARGLIST_H
#define INC_ARGLIST_H
/// Class: ArgList
/// Hold a list of string arguments, can be set from a delimited list
/// with SetList or arguments can be added individually with Add. 
/// Arguments can be accessed with the various getX routines,
/// where X is specific to certain types, e.g. getNextDouble returns
/// the next double, getNextMask returns an atom mask expression (i.e.
/// it has :, @, % characters etc).
class ArgList {
    char **arglist;
    char *marked;
    int nargs;
    char *argline;
    int debug;
  public:
    ArgList();
    ~ArgList();
    void SetDebug(int);

    int SetList(char*, const char*);
    ArgList *Copy();
    void Add(char *);

    void print();
    char *ArgLine();
    char *Arg(int);
    bool ArgIs(int,const char*);

    // In some cases (e.g. lines from input file) the first argument
    // is called the Command.
    char *Command();
    int CommandIs(const char *);
    int CommandIs(const char *,int);

    char *getNextString();
    void CheckForMoreArgs();
    char *getNextMask();
    int getNextInteger(int );
    double getNextDouble(double);
    char *getKeyString(const char *, char *);
    int getKeyIndex(char *);
    int getKeyInt(const char *, int);
    double getKeyDouble(const char *, double);
    int hasKey(const char *);

    int Contains(const char *);
    ArgList *SplitAt(const char *);
    ArgList *RemainingArgs();
    int ReplaceArg(int, char *);
    char *CopyArg(int);
    void Reset();
    void ResetAll();

    int Nargs() { return nargs; }
};
#endif
