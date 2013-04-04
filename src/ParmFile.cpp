#include "ParmFile.h"
#include "CpptrajStdio.h"
// All ParmIO classes go here
#include "Parm_Amber.h"
#include "Parm_PDB.h"
#include "Parm_Mol2.h"
#include "Parm_CharmmPsf.h"

const ParmFile::ParmToken ParmFile::ParmArray[] = {
  { AMBERPARM,    "amber", Parm_Amber::Alloc     },
  { PDBFILE,      "pdb",   Parm_PDB::Alloc       },
  { MOL2FILE,     "mol2",  Parm_Mol2::Alloc      },
  { CHARMMPSF,    "psf",   Parm_CharmmPsf::Alloc },
  { UNKNOWN_PARM, 0,       0                     }
};

// ParmFile::Read()
int ParmFile::Read(Topology& Top, std::string const& fname, bool bondsearch, int debugIn) 
{
  // Set up parm for reading
  CpptrajFile basicParm;
  int err = basicParm.SetupRead(fname, debugIn);
  if (err != 0) {
    mprinterr("Error: Could not set up parm file %s for reading.\n",fname.c_str());
    return 1;
  }
  
  // Loop over all parm formats
  for ( TokenPtr token = ParmArray; token->Alloc != 0; ++token ) {
    ParmIO* parmio = (ParmIO*)token->Alloc();
    parmio->SetDebug( debugIn );
    if ( parmio->ID_ParmFormat( basicParm ) ) {
      // Read this format
      err = parmio->ReadParm( basicParm.Filename().Full(), Top);
      parmName_ = basicParm.Filename();
      // Perform setup common to all parm files.
      if (err == 0) 
        err = Top.CommonSetup(bondsearch);
      else
        mprinterr("Error reading parm file %s\n", basicParm.Filename().full());
      delete parmio;
      return err;
    }
    delete parmio;
  }
  mprinterr("Error: Read: Format of parm [%s] not recognized.\n",fname.c_str());
  return 1;
}

// ParmFile::Write()
int ParmFile::Write(Topology const& Top, std::string const& fname, ParmFormatType fmt, int debugIn)
{
  // Loop over all parm formats
  for ( TokenPtr token = ParmArray; token->Alloc != 0; ++token ) {
    if ( token->Type == fmt ) {
      ParmIO* parmio = (ParmIO*)token->Alloc();
      parmio->SetDebug( debugIn );
      int err = parmio->WriteParm( fname, Top );
      if (err != 0 ) 
        mprinterr("Error writing parm file %s\n",fname.c_str());
      delete parmio;
      return err;
    }
  }
  mprinterr("Error: Write: Format of parm [%s] not recognized.\n",fname.c_str());
  return 1;
}
