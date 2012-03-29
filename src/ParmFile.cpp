#include <cstdio> // sscanf
#include <cstring> // strncmp, strlen
#include "ParmFile.h"
#include "PDBfileRoutines.h"
#include "Mol2FileRoutines.h"
#include "CpptrajStdio.h"
// All ParmIO classes go here
#include "Parm_Amber.h"
#include "Parm_PDB.h"
#include "Parm_Mol2.h"
#include "Parm_CharmmPsf.h"

// CONSTRUCTOR
ParmFile::ParmFile() {
  debug_ = 0;
}

// ParmFile::SetDebug() 
void ParmFile::SetDebug(int debugIn) {
  debug_ = debugIn;
}

// ParmFile::ID_ParmFormat()
ParmFile::ParmFormatType ParmFile::ID_ParmFormat(ParmIO &parm) {
  const int BUF_SIZE = 83;
  char buffer1[BUF_SIZE];
  char buffer2[BUF_SIZE];

  // Initialize buffer to NULL
  memset(buffer1,' ',BUF_SIZE);
  memset(buffer2,' ',BUF_SIZE);
  buffer1[0]='\0';
  buffer2[0]='\0';
  // Open Parm
  if (parm.OpenFile())
    return UNKNOWN_PARM;
  parm.IO->Gets(buffer1, BUF_SIZE);
  parm.IO->Gets(buffer2, BUF_SIZE);
  parm.CloseFile();

  // If both lines have PDB keywords, assume PDB
  if (isPDBkeyword(buffer1) && isPDBkeyword(buffer2))
  {
    if (debug_>0) mprintf("  PDB file\n");
    return PDBFILE;
  }

  // If either buffer contains a TRIPOS keyword assume Mol2
  // NOTE: This will fail on tripos files with extensive header comments.
  //       A more expensive check for mol2 files is below.
  if ( IsMol2Keyword(buffer1) || IsMol2Keyword(buffer2) )
  {
    if (debug_>0) mprintf("  TRIPOS MOL2 file\n");
    return MOL2FILE;
  }

  // If the %VERSION and %FLAG identifiers are present, assume amber parm
  if (strncmp(buffer1,"%VERSION",8)==0 && strncmp(buffer2,"%FLAG",5)==0) {
    if (debug_>0) mprintf("  AMBER TOPOLOGY file\n");
    return AMBERPARM;
  }

  // If the first 3 chars are P S F, assume charmm PSF
  if (strncmp(buffer1,"PSF",3)==0) {
    if (debug_>0) mprintf("  CHARMM PSF file\n");
    return CHARMMPSF;
  }

  // If first line is 81 bytes and the second line has 12 numbers in
  // 12I6 format, assume old-style Amber topology
  // NOTE: Could also be less than 81? Only look for 12 numbers?
  int iamber[12];
  if ((int)strlen(buffer1)==81+(int)parm.IsDos()) {
    if ( sscanf(buffer2,"%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i", iamber, iamber+1,
                iamber+2, iamber+3, iamber+4, iamber+5, iamber+6, 
                iamber+7, iamber+8, iamber+9, iamber+10, iamber+11) == 12 )
    {
      if (debug_>0) mprintf("  AMBER TOPOLOGY, OLD FORMAT\n");
      return OLDAMBERPARM;
    }
  }

  // ---------- MORE EXPENSIVE CHECKS ----------
  // Reopen and scan for Tripos mol2 molecule section
  // 0 indicates section found.
  if (parm.OpenFile())
    return UNKNOWN_PARM;
  if (Mol2ScanTo(parm.IO, MOLECULE)==0) {
    if (debug_>0) mprintf("  TRIPOS MOL2 file\n");
    parm.CloseFile();
    return MOL2FILE;
  }
  parm.CloseFile();

  // Unidentified file
  mprintf("  Warning: %s: Unknown topology format.\n",parm.BaseName());
  return UNKNOWN_PARM;
}

// ParmFile::Read()
int ParmFile::Read(AmberParm &parmOut, char *parm_filename,
                   bool bondsearch, bool molsearch) 
{
  ParmIO *parmio = NULL;
  ParmIO basicParm;

  int err = basicParm.SetupRead(parm_filename, debug_);
  if (err == 1) {
    mprinterr("Error setting up parm file %s for read.\n",parm_filename);
    return 1;
  }

  // Set Parm name
  parmOut.SetParmName( basicParm.BaseName(), basicParm.Name() );

  // Determine parm format
  ParmFormatType parmFormat = ID_ParmFormat( basicParm );

  switch ( parmFormat ) {
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

  // Place the basic file in the parm IO class
  parmio->ParmIO::operator=( basicParm );

  parmio->SetDebug(debug_);

  err = parmio->ReadParm( parmOut );
  if (err != 0) {
    mprinterr("Error reading parm file %s\n",parm_filename);
    delete parmio;
    return 1;
  }

  // Perform setup common to all parm files.
  parmOut.CommonSetup(bondsearch,molsearch);

  delete parmio;
  return 0;
}

// ParmFile::Write()
int ParmFile::Write(AmberParm &parmIn, char *parm_filename, ParmFormatType fmtIn) 
{
  ParmIO *parmio = NULL;
  ParmIO basicParm;

  int err = basicParm.SetupWrite(parm_filename, debug_);
  if (err == 1) {
    mprinterr("Error setting up parm file %s for write.\n",parm_filename);
    return 1;
  }

  // Set Parm name - NOTE: Make separate function?
  //parmIn.SetParmName( basicParm.basefilename, basicParm.filename );

  switch (fmtIn) {
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

  // Place the basic file in the parm IO class
  parmio->ParmIO::operator=( basicParm );

  parmio->SetDebug(debug_);

  err = parmio->WriteParm( parmIn );
  if (err != 0) {
    mprinterr("Error writing parm file %s\n",parm_filename);
    delete parmio;
    return 1;
  }

  delete parmio;
  return 0;
}

