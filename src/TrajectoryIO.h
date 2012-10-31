#ifndef INC_TRAJECTORYIO_H
#define INC_TRAJECTORYIO_H
#include "Topology.h" // Box
#include "CpptrajFile.h"
#include "ArgList.h"
// Class: TrajectoryIO
/// Abstract base class for performing trajectory reading and writing.
class TrajectoryIO {
  public:
    virtual ~TrajectoryIO(); // virtual since this class is inherited.
    struct TrajInfo {
      int NreplicaDim;     ///< Number of replica dimensions if applicable.
      Box BoxInfo;         ///< Box info for first frame of trajectory.
      bool HasV;           ///< True if trajectory has velocity info.
      bool HasT;           ///< True if trajectory has temperature info.
      bool IsSeekable;     ///< True if trajectory is randomly seekable.
      std::string Title;   ///< Set to trajectory title.
    };
    typedef TrajectoryIO* (*AllocatorType)();
    // -----------===== Inherited functions =====-----------
    /// Return true if file format matches trajectory type.
    virtual bool ID_TrajFormat(CpptrajFile&) = 0; 
    /// Set up trajectory IO for READ/APPEND
    /** Called inside TrajectoryFile::SetupRead. The first argument is the 
      * Topology class that will be associated with this trajectory. The
      * remaining argument will be set to contain trajectory information:
      * box (if any), velocity, temperature, etc.
      * \return the number of frames in the underlying trajectory file. 
      * \return -1 if an error occurs.
      */
    virtual int setupTrajin(std::string const&, Topology *, TrajInfo&) = 0;
    /// Set up trajectory IO for WRITE 
    /** Called inside TrajectoryFile::WriteFrame on the first write call. Takes
      * as arguments the Topology class that will be associated with this 
      * trajectory, the expected number of frames to be written, and any
      * additional trajectory information. 
      */
    virtual int setupTrajout(std::string const&,Topology*,int,TrajInfo const&,bool) = 0; 
    /// Open traj, prepare for IO.
    virtual int openTraj() = 0; 
    /// Read a frame from trajectory
    /** Given a frame number, read that frame; return the coordinates in the 
      * first array, velocities in the second array, the box lengths/angles in 
      * the third array, and set the temperature in the last var.
      */
    virtual int readFrame(int,double*,double*,double*,double*) = 0;
    /// Read velocity information from a trajectory. 
    virtual int readVelocity(int, double*) = 0;
    /// Read replica indices from a trajectory. 
    virtual int readIndices(int,int*)  = 0;  
    /// Write a frame to trajectory
    /** Write to output trajectory. This routine is called from
      * TrajectoryFile::WriteFrame with the current action set
      * number, not the current output number, so it is up to
      * the TrajectoryIO object to keep track of what frame it is
      * writing. Vars are same as in readFrame.
      */
    virtual int writeFrame(int,double*,double*,double*,double) = 0;
    /// Close trajectory
    virtual void closeTraj() = 0; 
    /// Print information on what kind of trajectory this is.
    virtual void info() = 0; 
    /// Process arguments relevant to writing trajectory (optional)
    /** Process any arguments from the arg list that have to do with 
      * setting the trajectory up for writing. Called before setupTrajout, so 
      * none of the arguments should be parm-related. It is desireable that any 
      * changes made to the TrajectoryIO object from within this function are
      * implemented as functions that can be called independently if need be 
      * (e.g. setting the write mode for PDB files).
      */
    virtual int processWriteArgs(ArgList&) = 0; 
    /// Process arguments relevant to reading trajectory (optional)
    virtual int processReadArgs(ArgList&) = 0; 
}; 
#endif
