#ifndef INC_TRAJECTORYIO_H
#define INC_TRAJECTORYIO_H
/// Class: TrajectoryIO
/// TrajectoryIO will be the base class for performing trajectory reading
/// and writing that all formats will inherit. If the trajectory format
/// reads/writes from a file, a PtrajFile object should be passed in via
/// the SetFile function. 
/// The following functions can be implemented by the inheriting class:
///   setupRead(): Called inside TrajectoryFile::SetupRead. Takes as an 
///                argument the expected number of atoms in the trajectory. 
///                Returns the number of frames in the underlying trajectory 
///                file. Should set all variables (title, seekable, hasBox, 
///                boxAngle (only if hasBox), and hasTemperature.
///   setupWrite(): Called inside TrajectoryFile::WriteFrame on the first
///                 write call. Takes as an argument the expected number of
///                 atoms in the trajectory. If any additional parm info is
///                 required it can be set prior to this call using
///                 SetParmInfo. 
///   openTraj(): Prepare trajectory for read/write
///   closeTraj(): Finish trajectory
///   readFrame(): Given a frame number, read that frame; return the
///                coordinates in the first array, the box lengths/angles in
///                the second array, and set the temperature in the last var.
///   writeFrame(): Write to output trajectory. This routine is called from
///                 TrajectoryFile::WriteFrame with the current action set
///                 number, not the current output number, so it is up to
///                 the TrajectoryIO object to keep track of what frame it is
///                 writing.
///   info(): Print information on what kind of trajectory this is.
#include "PtrajFile.h" 
class TrajectoryIO {
  protected:
    // **** This is defined so that arrays from AmberParm can be passed in
    //      via SetParmInfo, required for things like PDB/Mol2 writes. MUST
    //      match the definition in AmberParm or memory errors will occur.
    //      NOTE: Make this definition an include? 
    typedef char NAME[6];

    PtrajFile *tfile;   // Base file.
    char *title;        // Trajectory title.
    int debug;          // Debug level
    // Certain trajectory formats such as PDB and mol2 print out some
    // parm information in addition to coordinates. These variables are
    // only used for those formats.
    NAME *trajAtomNames;
    NAME *trajResNames;
    int *trajAtomsPerMol;
    int *trajResNums;
    double *trajCharges;
    double *trajRadii;
  public:
    bool seekable;      // True if can seek to frames in this traj.
    bool hasBox;        // True if the trajectory has box information.
    double boxAngle[3]; // Hold alpha, beta and gamma angles of box if hasBox.
    bool hasTemperature;// True if trajectory has temperature information.

    TrajectoryIO();
    virtual ~TrajectoryIO(); // virtual since this class is inherited.

    // Inherited functions 
    virtual int setupRead(int) { return -1; }
    virtual int setupWrite(int) { return 1; }
    virtual int openTraj() { return 1; }
    virtual int readFrame(int,double*,double*,double*) { return 1; }
    virtual int writeFrame(int,double*,double*,double) { return 1; }
    virtual void closeTraj() { return; }
    virtual void info() { return; }
  
    void SetFile(PtrajFile *);
    void SetTitle(char *);
    bool FilenameIs(char *);
    void SetDebug(int);
    void SetParmInfo(NAME*, NAME*, int *, int *, double *, double *);

    FileFormat TrajFormat() {return tfile->fileFormat;}
}; 
#endif
