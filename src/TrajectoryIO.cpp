// TrajectoryIO
#include "TrajectoryIO.h"
#include <cstring>
#include <cstdlib>
#include "CpptrajStdio.h"

// CONSTRUCTOR
TrajectoryIO::TrajectoryIO() {
  debug = 0;
  tfile=NULL;
  title=NULL;

  seekable=false;
  hasBox=false;
  boxAngle[0]=0.0;
  boxAngle[1]=0.0;
  boxAngle[2]=0.0;
  hasTemperature=false;

  trajAtomNames=NULL;
  trajResNames=NULL;
  trajAtomsPerMol=NULL;
  trajCharges=NULL;
  trajRadii=NULL;
}

// DESTRUCTOR
TrajectoryIO::~TrajectoryIO() {
  if (tfile!=NULL) delete tfile;
  if (title!=NULL) free(title);
}

/* TrajectoryIO::SetDebug()
 * Set debug level.
 */
void TrajectoryIO::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0) mprintf("TrajectoryIO debug level set to %i\n",debug);
}

/* TrajectoryIO::SetFile()
 * Set the internal file pointer.
 */
void TrajectoryIO::SetFile(PtrajFile *tfileIn) {
  tfile = tfileIn;
}

/* TrajectoryIO::FilenameIs()
 * Return true if trajectory full path filename matches input.
 */
bool TrajectoryIO::FilenameIs(char *filenameIn) {
  if (filenameIn==NULL) {
    mprinterr("Error: CheckFilename: Called with NULL filename.\n");
    return false;
  }
  if (strcmp(tfile->filename, filenameIn)==0) return true;
  return false;
}

/* TrajectoryIO::SetTitle()
 * Set title for this trajfile.
 */
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

/* TrajectoryIO::SetParmInfo()
 * Certain trajectory formats like PDB and Mol2 require some parm info
 * during writes. Set here.
 */
void TrajectoryIO::SetParmInfo(NAME *names, NAME *resnames, int *atomsPerMol, 
                               int *resnums, double *charge, double *radii) {
  //mprintf("SETTING PARM INFO FOR TRAJIO\n"); // DEBUG
  trajAtomNames = names;
  trajResNames = resnames;
  trajAtomsPerMol = atomsPerMol;
  trajResNums = resnums;
  trajCharges = charge;
  trajRadii = radii;
  // DEBUG
  //mprintf("Atomname 0 %s\n",trajAtomNames[0]);
  //mprintf("Resname 0  %s\n",trajResNames[0]);
  //mprintf("AtPerMol 0 %i\n",trajAtomsPerMol[0]);
  //mprintf("ResNums 0  %i\n",trajResNums[0]);
  //mprintf("Charge 0   %lf\n",trajCharges[0]);
}

