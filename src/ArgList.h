#ifndef INC_ARGLIST_H
#define INC_ARGLIST_H
// ArgList

class ArgList {
    char **arglist;
    char *marked;
    int nargs;
    char *argline;
  public:
    ArgList();
    ArgList(char *, const char *);
    ~ArgList();
    ArgList *Copy();
    void Add(char *);
    void print();
    char *ArgLine();
    char *Arg(int);

    char *Command();
    int CommandIs(const char *);
    int CommandIs(const char *,int);
    char *getNextString();
    void CheckForMoreArgs();
    char *getNextMask();
    int getNextInteger(int );
    char *getKeyString(const char *, char *);
    int getKeyIndex(char *);
    int getKeyInt(const char *, int);
    double getKeyDouble(const char *, double);
    int hasKey(const char *);
    ArgList *SplitAt(const char *);
    int ReplaceArg(int, char *);
    char *CopyArg(int);
    void Reset();
    void ResetAll();
};
#endif
