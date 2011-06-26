#ifndef INC_TRAJECTORYIO_H
#define INC_TRAJECTORYIO_H
/// Class: TrajectoryIO
/// TrajectoryIO will be the base class for performing trajectory reading
/// and writing that all formats will inherit. If the trajectory format
/// reads/writes from a file, a PtrajFile object should be passed in via
/// the SetFile function. 
/// The following functions should be implemented by the inheriting class:
///   setupRead(): Called inside TrajectoryFile::SetupRead. Takes as an 
///                argument the expected number of atoms in the trajectory. 
///                Returns the number of frames in the underlying trajectory 
///                file. Should set all variables (title, seekable, hasBox, 
///                boxAngle (only if hasBox), and hasTemperature.
///   openTraj(): Prepare trajectory for read/write
///   readFrame(): Given a frame number, read that frame; return the
///                coordinates in the first array, the box lengths/angles in
///                the second array, and set the temperature in the last var.
#include "PtrajFile.h" 
class TrajectoryIO {
  protected:
    PtrajFile *tfile;   // Base file.
    char *title;        // Trajectory title.
    int debug;          // Debug level
    // Certain trajectory formats such as PDB and mol2 print out some
    // parm information in addition to coordinates. These variables are
    // only used for those formats.
    char **trajAtomNames;
    char **trajResNames;
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
    virtual int setupRead(int) { return 1; }
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
    void SetParmInfo(char **, char **, int *, int *, double *, double *);

    FileFormat TrajFormat() {return tfile->fileFormat;}
}; 
#endif
