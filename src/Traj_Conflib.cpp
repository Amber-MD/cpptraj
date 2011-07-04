// Conflib
#include "Traj_Conflib.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Conflib::Conflib() { 
  energy=0.0;
  radGyr=0.0;
  timesFound=0;
  conflibAtom=0;
}

// DESTRUCTOR
Conflib::~Conflib() { 
  //fprintf(stderr,"Conflib Destructor.\n");
}
//------------------------------------------------------------------------
/* Conflib::closeTraj()
 */
void Conflib::closeTraj() {
  tfile->CloseFile();
}

/* Conflib::openTraj()
 */
int Conflib::openTraj() {
  if (tfile->OpenFile()) return 1;
  return 0;
}

/* Conflib::setupRead()
 */
int Conflib::setupRead(int natom) {
  long unsigned int confFrame;
  int Frames = 0;

  // Conflib is double,double,int,natom*3*double
  confFrame = (((natom * 3) + 2) * sizeof(double)) + sizeof(int);
  Frames = (int) (tfile->file_size / confFrame);

  if ( (tfile->file_size % confFrame) != 0 ) {
    mprintf("Warning: %s: Could not accurately predict # frames. This usually \n",
            tfile->filename);
    mprintf("         indicates a corrupted trajectory. Will attempt to read %i frames.\n",
            Frames);
  }
  conflibAtom = natom;
  return Frames;
}

/* Conflib::readFrame()
 */
int Conflib::readFrame(int set, double *X, double *V,double *box, double *T) {

  if (tfile->IO->Read(&energy,sizeof(double),1) < 0) return 1;
  tfile->IO->Read(&radGyr,sizeof(double),1);
  tfile->IO->Read(&timesFound,sizeof(int),1);
  tfile->IO->Read(X,sizeof(double),conflibAtom*3); 

  if (debug>0) mprinterr("CONFLIB %10i: %10.4lf %10.4lf %6i %10.4lf %10.4lf %10.4lf\n",
                         set, energy, radGyr, timesFound, X[0], X[1], X[2]);
  return 0;
}

/* Conflib::setupWrite()
 */
int Conflib::setupWrite(int natom ) {
  mprintf("Error: conflib writes not yet implemented.\n");
  return 1;
}

/* Conflib::setupWrite()
 */
int Conflib::writeFrame(int set, double *X, double *V,double *box, double T) {
  mprintf("Error: conflib writes not yet implemented.\n");
  return 1;
}

/* Conflib::info()
 */
void Conflib::info() {
  mprintf("is an LMOD conflib file");
}

