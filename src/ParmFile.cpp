#include "ParmFile.h"
#include "CpptrajStdio.h"
// All ParmIO classes go here
#include "Parm_Amber.h"
#include "Parm_PDB.h"
//#include "Parm_Mol2.h"
//#include "Parm_CharmmPsf.h"

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
    //case MOL2FILE    : parmio = new Mol2ParmFile(); break;
    //case CHARMMPSF   : parmio = new CharmmPsfParmFile(); break;
    default :
      mprinterr("Error: Parm file format not recognized (%i).\n",(int)parmFormat);
      return NULL;
  }
  return parmio;
}

// ParmFile::Read()
int ParmFile::Read(Topology &Top, char *fname, bool bondsearch, bool molsearch) {
  ParmIO *parmio = NULL;
  ParmIO basicParm;

  // Set up parm for reading
  int err = basicParm.SetupRead(fname, debug_);
  
  // Loop over all parm formats 
  for (int fmtidx = (int)UNKNOWN_PARM + 1; fmtidx != (int)NPARM; fmtidx++) {
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
          mprinterr("Error reading parm file %s\n",fname);
          delete parmio;
          return 1;
        }
        // Perform setup common to all parm files.
        Top.CommonSetup(bondsearch,molsearch);
        // Set base filename for later retrieval
        basename_.assign( parmio->BaseName() );
        delete parmio;
        return 0;
      }
      delete parmio;
    }
  }
    
  return 1;
}

// ParmFile::Write()
int ParmFile::Write( Topology &Top, std::string &name, ParmFormatType fmt) {
  return 1;
}

int ParmFile::Write( Topology &Top, char *nameIn, ParmFormatType fmt) {
  std::string name( nameIn );
  return Write(Top, name, fmt);
}

