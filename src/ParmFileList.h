#ifndef INC_PARMFILELIST_H
#define INC_PARMFILELIST_H
/// Class: ParmFileList
/// Holds a list of parameter files. Can either add new parm files
/// by filename, or add existing files by address. Search for parm
/// files in a list by index or full/base filename.
#include "AmberParm.h"
#include "ArgList.h"
class ParmFileList {
    std::vector<AmberParm*> ParmList;
    int Nparm;
    int debug;
    bool hasCopies;  // true: List contains addresses of parm files, do not delete
    bool bondsearch; // true: When parm is opened bond info will be filled in
    bool molsearch;  // true: When parm is opened molecule info will be filled in

  public:

    ParmFileList();
    ~ParmFileList();

    void SetDebug(int);
    int CheckCommand(ArgList *);
    int AddParmFile(char *);
    int AddParm(AmberParm *);
    AmberParm *GetParm(int);
    AmberParm *GetParm(ArgList &);
    int GetParmIndex(char *);
    int ReplaceParm(int, AmberParm *);
    void Print();
};
#endif
