#ifndef INC_COORDFILELIST_H
#define INC_COORDFILELIST_H
/// Class: CoordFileList
/// Base class for input trajectories (trajin), output trajectories (trajout),
/// and reference coordinates (reference). Each class that inherits this should
/// set fileAccess, which is the default access for files added to the list 
/// (e.g. for TrajinList, fileAccess is READ, etc). Each inheriting class will
/// also provide its own implementation of Add, which will Add files to the
/// list and set them up if appropriate.
#include <list>
#include "TrajectoryFile.h" // TrajectoryIO, AmberParm, ArgList, ProgressBar 
class CoordFileList {
  protected:
    std::list<TrajectoryFile*> trajList;
    std::list<TrajectoryFile*>::iterator currentTraj;
    AccessType fileAccess;      // READ/WRITE/APPEND, set in constructor/Add(write)
    int debug;                  // Debug level

    bool FilenameInUse(char *); // Return true if filename exists in list
  public:
    CoordFileList();
    virtual ~CoordFileList();   // Virtual since this class is inherited.

    void SetDebug(int);         // Set debug level
    void Info(int,int);         // Print information about all trajs in list

    // Trajectory List Functions
    void Begin();
    TrajectoryFile *NextTraj();

    // Inherited Functions
    virtual int Add(char *, ArgList *, AmberParm *) { return 1; }
};
#endif    
