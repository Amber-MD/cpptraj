// Parm_Pdb.cpp
#include <cstring> // memcpy, memset, strncpy, strcpy, strlen
#include "Parm_PDB.h"
#include "CpptrajStdio.h"
#include "PDBfileRoutines.h"

// PdbParmFile::ReadParm()
/** Open the PDB file specified by filename and set up topology data.
  * Mask selection requires natom, nres, names, resnames, resnums.
  */
int PdbParmFile::ReadParm(AmberParm &parmOut, CpptrajFile &parmfile) {
  char buffer[256];
  int atomInLastMol = 0;
  int lastResnum = -1;
  int currentResnum;

# ifdef USE_CHARBUFFER
  // TEST: Close and reopen buffered.
  parmfile.CloseFile();
  parmfile.OpenFileBuffered();
# endif

  mprintf("    Reading PDB file %s as topology file.\n",parmOut.parmName);
  memset(buffer,' ',256);
  // Initial readthrough to get # atoms, residues, molecules
# ifdef USE_CHARBUFFER
  while ( parmfile.Gets(buffer,256) == 0 )
# else
  while ( parmfile.IO->Gets(buffer,256)==0 )
# endif
  {
    // If ENDMDL or END is reached stop reading
    if ( strncmp(buffer,"END",3)==0) break;
    // If TER increment number of molecules
    if ( strncmp(buffer,"TER",3)==0) {
      ++parmOut.molecules;
      atomInLastMol = parmOut.natom;
    // Skip all other non-ATOM/HETATM keywords
    } else if (isPDBatomKeyword(buffer)) {
      // If this residue number is different than the last, increment nres
      currentResnum = pdb_resnum(buffer);
      if (currentResnum != lastResnum) {
        ++parmOut.nres;
        lastResnum = currentResnum;
        // If new residue and HETATM consider it a different molecule as well
        if (strncmp(buffer,"HETATM",6)==0) {
          // If HETATM immediately preceded by a TER card atomsPerMol has
          // just been set, so would be calling with 0. No need to call.
          if ( (parmOut.natom - atomInLastMol) != 0) {
            ++parmOut.molecules;
            atomInLastMol = parmOut.natom;
          }
        }
      }
      // Clear the buffer
      memset(buffer,' ',256);
      ++parmOut.natom;
    }
  }

  // If a TER card has been read and we are setting up the number of molecules,
  // finish up info on the last molecule read.
  if (parmOut.molecules>0 && (parmOut.natom - atomInLastMol) != 0)
    ++parmOut.molecules;

  // If no atoms, probably issue with PDB file
  if (parmOut.natom<=0) {
    mprintf("Error: No atoms in PDB file.\n");
    return 1;
  }
  if (debug>0) 
    mprintf("\tPDB contains %i atoms, %i residues, %i molecules.\n",
            parmOut.natom,parmOut.nres,parmOut.molecules);

  // Allocate memory
  parmOut.names = new NAME[ parmOut.natom ];
  parmOut.resnames = new NAME[ parmOut.nres ];
  parmOut.resnums = new int[ parmOut.nres ];
  parmOut.parmCoords = new double[ parmOut.natom * 3 ];
  if (parmOut.molecules > 0)
    parmOut.atomsPerMol = new int[ parmOut.molecules ];

  // Set / Reset indices 
  int iatm = 0;
  int ires = 0;
  int imol = 0;
  double *crdidx = parmOut.parmCoords;
  atomInLastMol = 0;
  lastResnum = -1;

  // Second readthrough to read in values
# ifdef USE_CHARBUFFER
  parmfile.Rewind();
  while ( parmfile.Gets(buffer,256) == 0 ) 
# else
  parmfile.IO->Rewind();
  while ( parmfile.IO->Gets(buffer,256)==0 ) 
# endif
  {
    // If ENDMDL or END is reached stop reading
    if ( strncmp(buffer,"END",3)==0) break;
    // If TER increment number of molecules and continue
    if ( strncmp(buffer,"TER",3)==0) {
      parmOut.atomsPerMol[imol++] = iatm - atomInLastMol;
      atomInLastMol = iatm;
    // Skip all other non-ATOM records
    } else if (isPDBatomKeyword(buffer)) {
      // Detect and remove trailing newline
      //bufferLen = strlen(buffer);
      //if (buffer[bufferLen-1] == '\n') buffer[bufferLen-1]='\0';
      // Copy atom name
      // Leading whitespace will automatically be trimmed.
      // Asterisks are replaced with single quote
      pdb_name(buffer, (char*)parmOut.names[iatm]);
      // Copy atom coords
      pdb_xyz(buffer,crdidx);
      crdidx+=3;

      // If this residue number is different than the last, copy 
      // residue name and atom number.
      currentResnum = pdb_resnum(buffer);
      if (currentResnum != lastResnum) {
        // Leading whitespace will automatically be trimmed.
        // Asterisks will be replaced with prime char
        pdb_resname(buffer, (char*)parmOut.resnames[ires]);
        if (debug>3) mprintf("\tPDBRes %i [%s]\n",ires+1,parmOut.resnames[ires]);
        parmOut.resnums[ires] = iatm;
        ++ires;
        lastResnum = currentResnum;
        // If new residue and HETATM consider it a different molecule as well
        if (strncmp(buffer,"HETATM",6)==0) {
          // If HETATM immediately preceded by a TER card atomsPerMol has
          // just been set, so would be calling with 0. No need to call.
          if ( (iatm - atomInLastMol) != 0) {
            parmOut.atomsPerMol[imol++] = iatm - atomInLastMol;
            atomInLastMol = iatm;
          }
        }
  
      // If residue number hasnt changed check if the current atom name has 
      // already been used in this residue.
      // NOTE: At this point ires has been incremented. Want ires-1.
      //       iatm is the current atom.
      } else {
        for (int atom = parmOut.resnums[ires-1]; atom < iatm; atom++) {
          if ( strcmp(parmOut.names[iatm], parmOut.names[atom])==0 ) {
            mprintf("\tWarning: Duplicate atom name in residue [%s]:%i [%s]:%i\n",
                    parmOut.resnames[ires-1],ires,parmOut.names[iatm],iatm+1);
          }
        }
      }
      // Clear the buffer
      memset(buffer,' ',256);
      ++iatm;
    } // END if atom/hetatm keyword
  } // END second read in parmfile

  // If a TER card has been read and we are setting up the number of molecules,
  // finish up info on the last molecule read.
  if (parmOut.molecules>0) {
    if (iatm - atomInLastMol != 0)
      parmOut.atomsPerMol[imol++] = iatm - atomInLastMol;
    // DEBUG
    if (debug>1) {
      mprintf("\tPDB: Atoms Per Molecule:\n");
      for (int mol=0; mol < parmOut.molecules; mol++) 
        mprintf("\t     %8i %8i\n",mol+1,parmOut.atomsPerMol[mol]);
    }
  }

  // No box for PDB - maybe change later to include unit cell info?
  parmOut.boxType = NOBOX;
  
  return 0;
}

