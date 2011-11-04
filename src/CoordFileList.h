#ifndef INC_COORDFILELIST_H
#define INC_COORDFILELIST_H
#include <list>
#include "TrajectoryFile.h" // TrajectoryIO, AmberParm, ArgList, ProgressBar 
// Class: CoordFileList
/// Base class for trajectory lists.
/** Used as base for input trajectories (trajin), output trajectories (trajout),
  * and reference coordinates (reference). Each class that inherits this should
  * set fileAccess, which is the default access for files added to the list 
  * (e.g. for TrajinList, fileAccess is READ, etc). Each inheriting class will
  * also provide its own implementation of Add, which will Add files to the
  * list and set them up if appropriate.
  */
class CoordFileList {
  protected:
    /// List of trajectories
    std::list<TrajectoryFile*> trajList;
    /// Iterator pointing to the current trajectory
    std::list<TrajectoryFile*>::iterator currentTraj;
    AccessType fileAccess;      ///< READ/WRITE/APPEND, set in constructor/Add(write)
    int debug;                  ///< Debug level

    bool FilenameInUse(char *); ///< Return true if filename exists in list
  public:
    CoordFileList();
    virtual ~CoordFileList();   // Virtual since this class is inherited.

    void SetDebug(int);         ///< Set debug level
    void Info(int,int);         ///< Print information about all trajs in list

    // Trajectory List Functions
    /// Position iterator at first trajectory
    void Begin();
    /// Return the current trajectory in the list, increment iterator
    TrajectoryFile *NextTraj();

    // Inherited Functions
    /// Add a trajectory to the list
    virtual int Add(char *, ArgList *, AmberParm *) { return 1; }
};
#endif    
