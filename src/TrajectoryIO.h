#ifndef INC_TRAJECTORYIO_H
#define INC_TRAJECTORYIO_H
#include "AmberParm.h" // BoxType
#include "CpptrajFile.h"
#include "ArgList.h"
// Class: TrajectoryIO
/// Base class for performing trajectory reading and writing.
/** This is the base class that all formats will inherit. If the trajectory format
  * reads/writes from a file, a CpptrajFile object should be passed in via
  * the SetFile function. 
  * The following functions can be implemented by the inheriting class:
  * setupRead, setupWrite, openTraj, closeTraj, readFrame, writeFrame, 
  * info, and processWriteArgs (optional).
  */
class TrajectoryIO {
  protected:
    CpptrajFile *tfile; ///< Base file.
    char *title;        ///< Trajectory title.
    int debug;          ///< Debug level
  public:
    bool seekable;      ///< True if can seek to frames in this traj.
    bool hasBox;        ///< True if the trajectory has box information.
    double boxAngle[3]; ///< Hold alpha, beta and gamma angles of box if hasBox.
    bool hasTemperature;///< True if trajectory has temperature information.
    bool hasVelocity;   ///< True if trajectory has velocity information.

    TrajectoryIO();
    virtual ~TrajectoryIO(); // virtual since this class is inherited.

    // -----===== Inherited functions =====-----
    /// Set up trajectory IO for READ/APPEND
    /** Called inside TrajectoryFile::SetupRead. Takes as an argument the 
      * AmberParm class that will be associated with this trajectory. Returns 
      * the number of frames in the underlying trajectory file. Should set all 
      * variables (title, seekable, hasBox, boxAngle (only if hasBox), 
      * hasTemperature, and hasVelocity. If an error occurs should return -1.
      */
    virtual int setupRead(AmberParm *) { return -1; }
    /// Set up trajectory IO for WRITE 
    /** Called inside TrajectoryFile::WriteFrame on the first write call. Takes
      * as an argument the AmberParm class that will be associated with this 
      * trajectory. 
      */
    virtual int setupWrite(AmberParm *) { return 1; }
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
      * setting the trajectory up for writing. Called before setupWrite, so 
      * none of the arguments should be parm-related. It is desireable that any 
      * changes made to the TrajectoryIO object from within this function are
      * implemented as functions that can be called independently if need be 
      * (e.g. setting the write  mode for PDB files).
      */
    virtual int processWriteArgs(ArgList *) { return 0; }
    /// Use the given CpptrajFile 
    void SetFile(CpptrajFile *);
    /// Set up a CpptrajFile
    int NewFile(char *, AccessType, FileFormat, FileType);
    /// Set the trajectory title
    void SetTitle(char *);
    /// Return true if trajectory filename is input
    bool FilenameIs(char *);
    /// Set debug level
    void SetDebug(int);
    /// For writes, indicate temperature info should be written if supported
    void SetTemperature();
    /// Return the trajectory format
    FileFormat TrajFormat();
    /// Return the filename
    char *Filename();
}; 
#endif
