#ifndef INC_TRAJECTORYIO_H
#define INC_TRAJECTORYIO_H
#include "Topology.h" // BoxType
#include "CpptrajFile.h"
#include "ArgList.h"
// Class: TrajectoryIO
/// Base class for performing trajectory reading and writing.
/** This is the base class that all formats will inherit. It inherits from 
  * CpptrajFile and so will be able to set up a FileIO class (IO) appropriate
  * to the underlying file type (e.g. standard, gzip, etc).  
  * The following functions can be implemented by the inheriting class:
  * setupTrajin, setupTrajout, openTraj, closeTraj, readFrame, writeFrame, 
  * info, and processWriteArgs (optional).
  */
class TrajectoryIO : protected CpptrajFile {
  public:
    TrajectoryIO();
    virtual ~TrajectoryIO(); // virtual since this class is inherited.
    TrajectoryIO(const TrajectoryIO&);
    TrajectoryIO& operator=(const TrajectoryIO&);
    TrajectoryIO& operator=(const CpptrajFile&);

    // -----------===== Inherited functions =====-----------
    /// Return true if file format matches trajectory type.
    virtual bool ID_TrajFormat() { return false; }
    /// Set up trajectory IO for READ/APPEND
    /** Called inside TrajectoryFile::SetupRead. Takes as an argument the 
      * Topology class that will be associated with this trajectory. Returns 
      * the number of frames in the underlying trajectory file. Should set all 
      * variables (title, seekable, hasBox, boxAngle (only if hasBox), 
      * hasTemperature, and hasVelocity. If an error occurs should return -1.
      */
    virtual int setupTrajin(Topology *) { return -1; }
    /// Set up trajectory IO for WRITE 
    /** Called inside TrajectoryFile::WriteFrame on the first write call. Takes
      * as arguments the Topology class that will be associated with this 
      * trajectory and the expected number of frames to be written. 
      */
    virtual int setupTrajout(Topology *, int) { return 1; }
    /// Open traj, prepare for IO.
    virtual int openTraj() { return 1; }
    /// Read a frame from trajectory
    /** Given a frame number, read that frame; return the coordinates in the 
      * first array, velocities in the second array, the box lengths/angles in 
      * the third array, and set the temperature in the last var.
      */
    virtual int readFrame(int,double*,double*,double*,double*) { return 1; }
    /// Write a frame to trajectory
    /** Write to output trajectory. This routine is called from
      * TrajectoryFile::WriteFrame with the current action set
      * number, not the current output number, so it is up to
      * the TrajectoryIO object to keep track of what frame it is
      * writing. Vars are same as in readFrame.
      */
    virtual int writeFrame(int,double*,double*,double*,double) { return 1; }
    /// Close trajectory
    virtual void closeTraj() { return; }
    /// Print information on what kind of trajectory this is.
    virtual void info() { return; }
    /// Process arguments relevant to writing trajectory (optional)
    /** Process any arguments from the arg list that have to do with 
      * setting the trajectory up for writing. Called before setupTrajout, so 
      * none of the arguments should be parm-related. It is desireable that any 
      * changes made to the TrajectoryIO object from within this function are
      * implemented as functions that can be called independently if need be 
      * (e.g. setting the write mode for PDB files).
      */
    virtual int processWriteArgs(ArgList *) { return 0; }
    // -----------------------------------------------------
    /// Set the trajectory title
    void SetTitle(char *);
    /// Set debug level - CpptrajFile has no SetDebug
    void SetDebug(int);
    /// For writes, indicate temperature info should be written if supported
    void SetTemperature();
    /// For writes, indicate box information should be written
    void SetBox();
    /// Check box info in traj against box info in parm
    int CheckBoxInfo(Topology*);

    /// Return true if the trajectory is seekable.
    bool Seekable()       { return seekable_;       }
    /// Return true if the trajectory has box coordinates.
    bool HasBox()         { return hasBox_;         }
    /// Return true if the trajectory has temperature info.
    bool HasTemperature() { return hasTemperature_; }
    /// Return true if the trajectory has velocity info.
    bool HasVelocity()    { return hasVelocity_;    }
  protected:
    std::string title_;  ///< Trajectory title.
    bool seekable_;      ///< True if can seek to frames in this traj.
    bool hasBox_;        ///< True if the trajectory has box information.
    double boxAngle_[3]; ///< Hold alpha, beta and gamma angles of box if hasBox.
    double boxLength_[3];///< Hold x, y and z lengths of box if hasBox.
    bool hasTemperature_;///< True if trajectory has temperature information.
    bool hasVelocity_;   ///< True if trajectory has velocity information.
}; 
#endif
