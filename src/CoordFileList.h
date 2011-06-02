#ifndef INC_COORDFILELIST_H
#define INC_COORDFILELIST_H
/*
 * CoordFileList
 * Base class for input trajectories (trajin), output trajectories (trajout),
 * and reference coordinates (reference). Each class that inherits this should
 * set Command, which is the keyword required to add trajectories to the list
 * (e.g. for TrajinList Command is "trajin" etc). Each inheriting class will
 * also provide its own implementation of Add, which will call ProcessArgList
 * to get arguments common to all traj list types. 
 */
#include <list>
#include "TrajFile.h" // ArgList.h

class CoordFileList : public std::list<TrajFile *> {
  protected:
    //char *trajfilename;               // Traj filename, from Arg list
    //AmberParm *P;                     // Parm, determined from Args (parm/parmindex)
    AccessType fileAccess;            // READ/WRITE/APPEND, set in constructor/Add(write)
    int debug;                        // Debug level
    std::list<TrajFile *>::iterator it; // Iterator for the list

    int CheckFilename(char *);
    //int ProcessArgList(ArgList *, ParmFileList *);

  public:
    CoordFileList();
    virtual ~CoordFileList();

    void SetDebug(int);
    TrajFile *SetupTrajectory(char *, AccessType, FileFormat, FileType);
    void Info(int);

    virtual int Add(ArgList *, AmberParm *) { return 1; }
};
#endif    
