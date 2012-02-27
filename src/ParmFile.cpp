#include "ParmFile.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"
// All ParmIO classes go here
#include "Parm_Amber.h"
#include "Parm_PDB.h"
#include "Parm_Mol2.h"
#include "Parm_CharmmPsf.h"

// CONSTRUCTOR
ParmFile::ParmFile() {
  debug = 0;
}

// ParmFile::SetDebug() 
void ParmFile::SetDebug(int debugIn) {
  debug = debugIn;
}

// ParmFile::Read()
int ParmFile::Read(AmberParm &parmOut, char *parm_filename,
                   bool bondsearch, bool molsearch) 
{
  CpptrajFile basicParm;
  ParmIO *parmio = NULL;

  int err = basicParm.SetupFile(parm_filename, READ, UNKNOWN_FORMAT, UNKNOWN_TYPE, debug);
  if (err == 1) {
    mprinterr("Error setting up parm file %s for read.\n",parm_filename);
    return 1;
  }

  // Set Parm name
  parmOut.SetParmName( basicParm.basefilename, basicParm.filename );

  switch (basicParm.fileFormat) {
    case OLDAMBERPARM: 
    case AMBERPARM   : parmio = new AmberParmFile(); break;
    case PDBFILE     : parmio = new PdbParmFile(); break;
    case MOL2FILE    : parmio = new Mol2ParmFile(); break;
    case CHARMMPSF   : parmio = new CharmmPsfParmFile(); break;
    default :
      mprinterr("Error: Parm file %s format not recognized.\n",parm_filename);
      return 1;
  }

  if (parmio==NULL) return 1;

  parmio->SetDebug(debug);

  basicParm.OpenFile();
  err = parmio->ReadParm( parmOut, basicParm );
  if (err != 0) {
    mprinterr("Error reading parm file %s\n",parm_filename);
    delete parmio;
    return 1;
  }
  basicParm.CloseFile();

  // Perform setup common to all parm files.
  parmOut.CommonSetup(bondsearch,molsearch);

  delete parmio;
  return 0;
}

// ParmFile::Write()
int ParmFile::Write(AmberParm &parmIn, char *parm_filename, FileFormat fmtIn) 
{
  CpptrajFile basicParm;
  ParmIO *parmio = NULL;

  int err = basicParm.SetupFile(parm_filename, WRITE, fmtIn, UNKNOWN_TYPE, debug);
  if (err == 1) {
    mprinterr("Error setting up parm file %s for write.\n",parm_filename);
    return 1;
  }

  // Set Parm name - NOTE: Make separate function?
  //parmIn.SetParmName( basicParm.basefilename, basicParm.filename );

  switch (basicParm.fileFormat) {
    //case OLDAMBERPARM: parmio = new OldAmberParmFile(); break;
    case AMBERPARM   : parmio = new AmberParmFile(); break;
    //case PDBFILE     : parmio = new PdbParmFile(); break;
    //case MOL2FILE    : parmio = new Mol2ParmFile(); break;
    //case CHARMMPSF   : parmio = new CharmmPsfParmFile(); break;
    default :
      mprinterr("Error: Parm file %s format not recognized.\n",parm_filename);
      return 1;
  }

  if (parmio==NULL) return 1;

  parmio->SetDebug(debug);

  basicParm.OpenFile();
  err = parmio->WriteParm( parmIn, basicParm );
  if (err != 0) {
    mprinterr("Error writing parm file %s\n",parm_filename);
    delete parmio;
    return 1;
  }
  basicParm.CloseFile();

  delete parmio;
  return 0;
}

