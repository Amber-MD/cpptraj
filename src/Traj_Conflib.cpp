// Conflib
#include "Conflib.h"
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
/* 
 * Conflib::close()
 */
void Conflib::close() {
  File->CloseFile();
}

/* 
 * Conflib::open()
 */
int Conflib::open() {

  if (File->OpenFile()) return 1;

  return 0;
}


/* 
 * Conflib::SetupRead()
 */
int Conflib::SetupRead() {
  long unsigned int confFrame;

  // Conflib is double,double,int,natom*3*double
  confFrame = (((P->natom * 3) + 2) * sizeof(double)) + sizeof(int);

  if ( (File->file_size % confFrame) == 0 ) {
    Frames = (int) (File->file_size / confFrame);
    stop = Frames;
  } else {
    mprintf("Warning: Conflib::SetupRead(): Could not predict # frames\n");
    mprintf("         Ensure that associated parm has correct # atoms.\n");
    mprintf("         File size=%lu confFrame=%lu\n",File->file_size,
            confFrame);
    Frames=-1;
    stop=-1;
    conflibAtom = P->natom;
    return 1;
  }

  return 0;
}

/* 
 * Conflib::getFrame()
 * Use conflibAtom instead of P->natom in case of stripped parmtop
 */
int Conflib::getFrame(int set) {

  if (File->IO->Read(&energy,sizeof(double),1) < 0) return 1;
  File->IO->Read(&radGyr,sizeof(double),1);
  File->IO->Read(&timesFound,sizeof(int),1);
  File->IO->Read(F->X,sizeof(double),conflibAtom*3); 

  return 0;
}

// Set up trajectory for either write or append
int Conflib::SetupWrite( ) {
  mprintf("Error: conflib writes not yet implemented.\n");
  return 1;
}

// Write a frame
int Conflib::writeFrame(int set) {
  mprintf("Error: conflib writes not yet implemented.\n");
  return 1;
}

/*
 * Info()
 */
void Conflib::Info() {
  mprintf("is an LMOD conflib file");
}
