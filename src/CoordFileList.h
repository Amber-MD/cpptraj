#ifndef INC_COORDFILELIST_H
#define INC_COORDFILELIST_H
/// Class: CoordFileList
/// Base class for input trajectories (trajin), output trajectories (trajout),
/// and reference coordinates (reference). Each class that inherits this should
/// set Command, which is the keyword required to add trajectories to the list
/// (e.g. for TrajinList Command is "trajin" etc). Each inheriting class will
/// also provide its own implementation of Add, which will call ProcessArgList
/// to get arguments common to all traj list types. 
#include <list>
#include "TrajectoryFile.h" // TrajectoryIO, AmberParm, ArgList, ProgressBar 
class CoordFileList : public std::list<TrajectoryFile*> {
  protected:
    AccessType fileAccess;      // READ/WRITE/APPEND, set in constructor/Add(write)
    int debug;                  // Debug level

    bool FilenameInUse(char *); // Return true if filename exists in list
  public:
    CoordFileList();
    virtual ~CoordFileList();   // Virtual since this class is inherited.

    void SetDebug(int);         // Set debug level
    void Info(int);             // Print information about all trajs in list

    // Inherited Functions
    virtual int Add(char *, ArgList *, AmberParm *) { return 1; }
};
#endif    
