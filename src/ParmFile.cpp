#include "ParmFile.h"
#include "CpptrajStdio.h"
// All ParmIO classes go here
#include "Parm_Amber.h"
#include "Parm_PDB.h"
#include "Parm_Mol2.h"
#include "Parm_CharmmPsf.h"

// CONSTRUCTOR
ParmFile::ParmFile() :
  debug_(0)
{}

// ParmFile::SetDebug() 
void ParmFile::SetDebug(int debugIn) {
  debug_ = debugIn;
}

/// Set up ParmIO class for the given format
ParmIO *ParmFile::SetupParmIO(ParmFormatType parmFormat) {
  ParmIO *parmio = NULL;
  switch ( parmFormat ) {
    case AMBERPARM   : parmio = new Parm_Amber(); break;
    case PDBFILE     : parmio = new Parm_PDB(); break;
    case MOL2FILE    : parmio = new Parm_Mol2(); break;
    case CHARMMPSF   : parmio = new Parm_CharmmPsf(); break;
    default :
      mprinterr("Error: Parm file format not recognized (%i).\n",(int)parmFormat);
      return NULL;
  }
  return parmio;
}

// ParmFile::Read()
int ParmFile::Read(Topology &Top, std::string const& fname, bool bondsearch) 
{
  ParmIO *parmio = NULL;
  CpptrajFile basicParm;

  // Set up parm for reading
  int err = basicParm.SetupRead(fname, debug_);
  if (err!=0) {
    mprinterr("Error: Could not set up parm file %s for reading.\n",fname.c_str());
    return 1;
  }
  
  // Loop over all parm formats 
  for (int fmtidx = (int)UNKNOWN_PARM + 1; fmtidx != (int)NPARM; fmtidx++) {
    //mprintf("DEBUG: Checking Format %i\n",fmtidx);
    // Set up parmio for this format
    parmio = SetupParmIO( (ParmFormatType)fmtidx );
    if (parmio!=NULL) {
      // Set parmio file information from basicParm
      parmio->ParmIO::operator=( basicParm );
      //parmio->SetDebug( debug_ ); // NOTE: Debug should be set from above assignment
      // Check if this file matches current parmio format
      if (parmio->ID_ParmFormat()) {
        // Read this format
        err = parmio->ReadParm( Top );
        if (err!=0) {
          mprinterr("Error reading parm file %s\n",fname.c_str());
          delete parmio;
          return 1;
        }
        // Perform setup common to all parm files.
        Top.CommonSetup(bondsearch);
        // Set base filename for later retrieval
        basename_ = basicParm.BaseFileName();
        delete parmio;
        return 0;
      }
      delete parmio;
    }
  }
  mprinterr("Error: Format of parm [%s] not recognized.\n",fname.c_str());
    
  return 1;
}

// ParmFile::Write()
int ParmFile::Write( Topology &Top, std::string const& name, ParmFormatType fmt) {
  ParmIO *parmio = NULL;
  CpptrajFile basicParm;
  // Set up basic parm file for write
  int err = basicParm.SetupWrite(name, debug_);
  if (err == 1) {
    mprinterr("Error setting up parm file %s for write.\n",name.c_str());
    return 1;
  }
  // Setup ParmIO for this format
  parmio = SetupParmIO( fmt );
  if (parmio == NULL) return 1;
  // Place the basic file in the ParmIO class
  // NOTE: Should also set debug level.
  parmio->ParmIO::operator=( basicParm );
  // TODO: Set basename?

  err = parmio->WriteParm( Top );
  if (err != 0) {
    mprinterr("Error writing parm file %s\n",name.c_str());
    delete parmio;
    return 1;
  }

  delete parmio;

  return 0;
}

