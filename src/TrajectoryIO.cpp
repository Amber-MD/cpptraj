// TrajectoryIO
#include "TrajectoryIO.h"
#include <cstring>
#include <cstdlib>
#include "CpptrajStdio.h"

// CONSTRUCTOR
TrajectoryIO::TrajectoryIO() {
  tfile=NULL;
  title=NULL;
  debug = 0;

  seekable=false;
  hasBox=false;
  boxAngle[0]=0.0;
  boxAngle[1]=0.0;
  boxAngle[2]=0.0;
  hasTemperature=false;
  hasVelocity=false;
}

// DESTRUCTOR
TrajectoryIO::~TrajectoryIO() {
  if (tfile!=NULL) delete tfile;
  if (title!=NULL) free(title);
}

// TrajectoryIO::SetDebug()
/** Set debug level. */
void TrajectoryIO::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0) mprintf("TrajectoryIO debug level set to %i\n",debug);
}

// TrajectoryIO::SetFile()
/** Set the internal file pointer.  */
void TrajectoryIO::SetFile(CpptrajFile *tfileIn) {
  tfile = tfileIn;
}

// TrajectoryIO::NewFile()
/** For use when TrajectoryIO has already been set up to a specific format,
  * set up the IO file. the given format MUST match the type for which
  * this file has been set up, otherwise bad things could happen.
  */
int TrajectoryIO::NewFile(char *filenameIn, AccessType accIn, FileFormat fmtIn, 
                          FileType typeIn) 
{
  tfile = new CpptrajFile();
  if (tfile->SetupFile(filenameIn,accIn,fmtIn,typeIn,debug)) {
    //mprinterr("    Error: Could not set up file %s.\n",tname);
    delete tfile;
    return 1;
  }
  return 0;
}

// TrajectoryIO::FilenameIs()
/** Return true if trajectory full path filename matches input. */
bool TrajectoryIO::FilenameIs(char *filenameIn) {
  if (filenameIn==NULL) {
    mprinterr("Error: CheckFilename: Called with NULL filename.\n");
    return false;
  }
  if (strcmp(tfile->filename, filenameIn)==0) return true;
  return false;
}

// TrajectoryIO::FileFormat()
/** Return the format of associated file. */
FileFormat TrajectoryIO::TrajFormat() {
  if (tfile==NULL) return UNKNOWN_FORMAT;
  return (tfile->fileFormat);
}

// TrajectoryIO::SetTitle()
/** Set title for this trajfile. */
void TrajectoryIO::SetTitle(char *titleIn) {
  size_t titleSize;

  //mprintf("DEBUG: Attempting to set title for %s: [%s]\n",trajfilename,titleIn);
  if (titleIn==NULL) return;
  titleSize = strlen(titleIn);
  //mprintf("       Title size is %i\n",titleSize);
  if (titleSize==0) {
    mprintf("Warning: SetTitle(): Title for %s is 0 length.\n",tfile->filename);
    return;
  }
  this->title = (char*) realloc(this->title, (titleSize+1) * sizeof(char));
  if (this->title==NULL) {
    mprinterr("Error: SetTitle(): Could not allocate memory for title of %s.\n",tfile->filename);
    return;
  }
  strcpy(this->title, titleIn);
  // Remove terminal newline if one exists
  if (this->title[titleSize-1]=='\n') this->title[titleSize-1]='\0';
}

// TrajectoryIO::SetTemperature()
/** Indicate temperature information should be read/written. */
void TrajectoryIO::SetTemperature() {
  hasTemperature = true;
}

// TrajectoryIO::Filename()
char *TrajectoryIO::Filename() {
  if (tfile==NULL) return NULL;
  return tfile->filename;
}
