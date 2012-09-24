#include "Trajin_Single.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Trajin_Single::Trajin_Single() :
  trajio_(0)
{}

// DESTRUCTOR
Trajin_Single::~Trajin_Single() {
  if (trajio_!=0) delete trajio_;
}

// Trajin_Single::SetupTrajRead()
int Trajin_Single::SetupTrajRead(std::string const& tnameIn, ArgList *argIn, Topology *tparmIn) 
{
  // Require a filename
  if (tnameIn.empty()) {
    mprinterr("Internal Error: Trajin_Single: No filename given.\n");
    return 1;
  }
  // Check and set associated parm file
  if ( SetTrajParm( tparmIn ) ) return 1;
  // Check that file exists
  if (!fileExists(tnameIn.c_str())) {
    mprinterr("Error: File %s does not exist.\n",tnameIn.c_str());
    return 1;
  }
  // Set up CpptrajFile
  CpptrajFile baseFile;
  if (baseFile.SetupRead( tnameIn, debug_ )) return 1;
  // Set file name and base file name
  SetFileNames( tnameIn, baseFile.BaseFileName() );
  // Detect format
  if ( (trajio_ = DetectFormat( baseFile )) == 0 ) {
    mprinterr("Error: Could not set up file %s for reading.\n", tnameIn.c_str());
    return 1;
  }
  // Set up the format for reading and get the number of frames.
  if (StartStopOffset( trajio_, argIn )) return 1;
  // Check how many frames will actually be read
  if (setupFrameInfo() == 0) return 1;
  return 0;
}

