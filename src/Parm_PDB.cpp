// Parm_Pdb.cpp
#include <cstring> // memcpy, memset, strncpy, strcpy, strlen
#include "Parm_PDB.h"
#include "CpptrajStdio.h"
#include "PDBfileRoutines.h"

// PdbParmFile::SetAtomsPerMolPDB()
/** Use in ReadParmPDB only, when TER is encountered or end of PDB file
  * update the atomsPerMol array. Take number of atoms in the molecule
  * (calcd as current #atoms - #atoms in previous molecule) as input. 
  * Check if the last residue is solvent; if so, set up solvent information.
  * \return the current number of atoms.
  */
int PdbParmFile::SetAtomsPerMolPDB(AmberParm &parmOut, int atomInLastMol) {
  int numAtoms = parmOut.natom - atomInLastMol;
  //mprintf("DEBUG:\tCalling SetAtomsPerMolPDB with %i\n",numAtoms);
  if (numAtoms<1) return 0;
  // Check if the current residue is a solvent molecule
  //mprintf("DEBUG: Checking if %s is solvent.\n",resnames[nres-1]);
  //if (nres>0 && IsSolventResname(resnames[nres-1])) {
  //  if (firstSolvMol==-1) {
  //    firstSolvMol = molecules + 1; // +1 to be consistent w/ Amber top
  //    finalSoluteRes = nres - 1;    // +1 to be consistent w/ Amber top
  //  }
  //}
  int *tempAPM = new int[ parmOut.molecules+1 ];
  memcpy(tempAPM, parmOut.atomsPerMol, parmOut.molecules * sizeof(int));
  delete[] parmOut.atomsPerMol;
  parmOut.atomsPerMol = tempAPM;
  parmOut.atomsPerMol[parmOut.molecules] = numAtoms;
  parmOut.molecules++;
  return parmOut.natom;
}


// PdbParmFile::ReadParm()
/** Open the PDB file specified by filename and set up topology data.
  * Mask selection requires natom, nres, names, resnames, resnums.
  */
int PdbParmFile::ReadParm(AmberParm &parmOut, CpptrajFile &parmfile) {
  char buffer[256];
  int bufferLen;  
  int currResnum;
  int atomInLastMol = 0;
  unsigned int crdidx = 0;

# ifdef USE_CHARBUFFER
  // TEST: Close and reopen buffered.
  parmfile.CloseFile();
  parmfile.OpenFileBuffered();
# endif

  mprintf("    Reading PDB file %s as topology file.\n",parmOut.parmName);
  currResnum=-1;
  memset(buffer,' ',256);
# ifdef USE_CHARBUFFER
  while ( parmfile.Gets(buffer,256) == 0 ) 
# else
  while ( parmfile.IO->Gets(buffer,256)==0 ) 
# endif
  {
    // If ENDMDL or END is reached stop reading
    if ( strncmp(buffer,"END",3)==0) break;
    // If TER increment number of molecules and continue
    if ( strncmp(buffer,"TER",3)==0) {
      atomInLastMol = SetAtomsPerMolPDB(parmOut, atomInLastMol);
      continue;
    }
    // Skip all other non-ATOM records
    if (isPDBatomKeyword(buffer)) {
      // Detect and remove trailing newline
      bufferLen = strlen(buffer);
      if (buffer[bufferLen-1] == '\n') buffer[bufferLen-1]='\0';

      // Allocate memory for atom name.
      NAME *tempname = new NAME[ parmOut.natom+1 ];
      memcpy(tempname, parmOut.names, parmOut.natom * sizeof(NAME));
      delete[] parmOut.names;
      parmOut.names = tempname;
      // Leading whitespace will automatically be trimmed.
      // Name will be wrapped if it starts with a digit.
      // Asterisks will be replaced with prime char
      pdb_name(buffer, (char*)parmOut.names[parmOut.natom]);

      // Allocate memory for coords
      double *tempcoord = new double[ (parmOut.natom+1)*3 ];
      memcpy(tempcoord, parmOut.parmCoords, (parmOut.natom*3)*sizeof(double) );
      delete[] parmOut.parmCoords;
      parmOut.parmCoords = tempcoord;
      pdb_xyz(buffer, parmOut.parmCoords + crdidx);
      crdidx+=3;
      // If this residue number is different than the last, allocate mem for new res
      if (currResnum!=pdb_resnum(buffer)) {
        NAME *temprname = new NAME[ parmOut.nres+1 ];
        memcpy(temprname, parmOut.resnames, parmOut.nres * sizeof(NAME));
        delete[] parmOut.resnames;
        parmOut.resnames = temprname;
        // Leading whitespace will automatically be trimmed.
        // Asterisks will be replaced with prime char
        pdb_resname(buffer, (char*)parmOut.resnames[parmOut.nres]);
        if (debug>3) mprintf("        PDBRes %i [%s]\n",parmOut.nres,parmOut.resnames[parmOut.nres]);
        int *temprnum = new int[ parmOut.nres+1 ];
        memcpy(temprnum, parmOut.resnums, parmOut.nres * sizeof(int));
        delete[] parmOut.resnums;
        parmOut.resnums = temprnum;
        parmOut.resnums[parmOut.nres]=parmOut.natom; 
        currResnum=pdb_resnum(buffer);
        parmOut.nres++;
        
        // If new residue and HETATM consider it a different molecule as well
        if (strncmp(buffer,"HETATM",6)==0) {
          // If HETATM immediately preceded by a TER card atomsPerMol has
          // just been set, so would be calling with 0. No need to call.
          if ( (parmOut.natom - atomInLastMol) != 0)
            atomInLastMol = SetAtomsPerMolPDB(parmOut, atomInLastMol);
        }
  
      // If residue number hasnt changed check for duplicate atom names in res
      // NOTE: At this point nres has been incremented. Want nres-1.
      //       natom is the current atom.
      } else {
        for (int atom=parmOut.resnums[parmOut.nres-1]; atom < parmOut.natom; atom++) {
          if ( strcmp(parmOut.names[parmOut.natom], parmOut.names[atom])==0 ) {
            mprintf("      Warning: Duplicate atom name in residue [%s]:%i [%s]:%i\n",
                    parmOut.resnames[parmOut.nres-1],parmOut.nres,parmOut.names[parmOut.natom],parmOut.natom+1);
          }
        }
      }
      // Clear the buffer
      memset(buffer,' ',256);

      parmOut.natom++;
    } // END if atom/hetatm keyword
  } // END read in parmfile

  // If a TER card has been read and we are setting up the number of molecules,
  // finish up info on the last molecule read.
  if (parmOut.molecules>0) {
    SetAtomsPerMolPDB(parmOut, atomInLastMol);
    // DEBUG
    if (debug>0) {
      //mprintf("\tPDB: firstSolvMol= %i\n",firstSolvMol);
      mprintf("\tPDB: finalSoluteRes= %i\n",parmOut.finalSoluteRes);
      if (debug>1) {
        mprintf("\tPDB: Atoms Per Molecule:\n");
        for (int atom=0; atom < parmOut.molecules; atom++) {
          mprintf("\t     %8i %8i\n",atom,parmOut.atomsPerMol[atom]);
        } 
      }
    }
  }

  // No box for PDB - maybe change later to include unit cell info?
  parmOut.boxType = NOBOX;

  if (debug>0) 
    mprintf("\tPDB contains %i atoms, %i residues, %i molecules.\n",
            parmOut.natom,parmOut.nres,parmOut.molecules);
  // If no atoms, probably issue with PDB file
  if (parmOut.natom<=0) {
    mprintf("Error: No atoms in PDB file.\n");
    return 1;
  }

  return 0;
}

