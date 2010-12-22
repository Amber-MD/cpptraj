#ifndef INC_TRAJFILE_H
#define INC_TRAJFILE_H
/* TrajFile is a base class that all other trajectory types inherit.
 * Once the trajectory subtype is determined its address will be assigned
 * to a TrajFile pointer.
 * e.g.: If the traj is an Amber Trajectory, TrajFile *T=&AmberTraj.
 * The trajectory can then be accessed through virtual members.
 */
#include <list> // For FrameRange
#include "Frame.h"
//#include "ArgList.h"
#include "AmberParm.h" // PtrajFile.h cstdio
//#include "PtrajFile.h" // cstdio

class TrajFile {
  protected:
    char *title;        // The trajectory title
    int outputStart;    // where output should begin
    int frameskip;      // Number of frames to skip while reading; if not seekable this is 1
    int targetSet;      // The next frame to be processed
    int seekable;       // =1 if the file can be randomly accessed
    int currentFrame;   // Current frame in trajectory
    int showProgress;   // If 1, show progressbar during traj processing
    int start;          // Frame to start processing
    int stop;           // Frame to end processing
    int offset;         // Number of frames to skip while processing

    // --== Inherited by child classes ==--
    virtual int open() { return 0; }        // Open the file, prepare for coord read/write
    virtual void close() {}                 // Close the file

  public:
    int debug;             // Level of debug information to print
    char *trajfilename;    // The base trajectory filename
    // NOTE: I hate that the following are public. Only necessary for REMD processing!!
    int Frames;            // Total number of frames in trajectory
    int total_read_frames; // Total number of frames that will be read
    int isBox;             // >0 means trajectory has box information

    std::list<int> *FrameRange; // list of frames to be written out
    int hasTemperature;    // 1 means trajectory has temperature information
    PtrajFile *File;       // Class that handles basic file IO
    AmberParm *P;          // Memory address of the associated parmfile
    Frame *F;              // Hold coordinates of the current frame
    int skip;              /* READ: If =1 do not process this input trajectory 
                              WRITE: If =1 this traj has been set up for write */

    TrajFile();            // Constructor
    virtual ~TrajFile();   // Destructor - virtual since this class is inherited.

    void SetTitle(char *);   // Set trajectory title.
    void PrintInfo(int);     // Print trajectory Information
    int setupFrameInfo(int); // Set actual start/stop based on total #frames and #threads 
    int Begin(int *, int);   /* Prepare traj for processing. Set output start value, calcd in 
                              * setupFrameInfo. Allocate memory for F. 
                              */
    int Begin();                 // Prepare trajectory for output
    int NextFrame(int*);         // Put the next target frame into F.
    void End();                  // Close trajectory and free F memory
    void progressBar();          // Display trajectory progress to screen
//    void progressBar2();         // Display trajectory progress to screen
   
    void SetArgs(int,int,int);   // Set the stop, start, and offset args from user input
    // --== Inherited by child classes ==--
    virtual int getFrame(int) { return 0; } // Read the next coord frame into F
    virtual int SetupRead()     { return 0; } // Set file up for reading
    virtual int SetupWrite()    { return 0; } // Set file up for writing
    virtual int writeFrame(int) { return 0; } // Write coords in frame F to file
    virtual void Info() { return; }           // Print information about this trajectory
};
#endif
