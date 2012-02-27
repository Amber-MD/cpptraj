// Parm_Mol2.cpp
#include "Parm_Mol2.h"
#include "Mol2FileRoutines.h"
#include "CpptrajStdio.h"
#include <cstdio> // sscanf
#include <cstring> // memcpy, strcpy

// Mol2ParmFile::ReadParm()
/** Read file as a Tripos Mol2 file. */
int Mol2ParmFile::ReadParm(AmberParm &parmOut, CpptrajFile &parmfile) {
  char buffer[MOL2BUFFERSIZE];
  int mol2bonds;
  int resnum, currentResnum;
  unsigned int crdidx = 0;
  char resName[5];

  currentResnum=-1;
  mprintf("    Reading Mol2 file %s as topology file.\n",parmOut.parmName);
  // Get @<TRIPOS>MOLECULE information
  if (Mol2ScanTo(&parmfile, MOLECULE)) return 1;
  //   Scan title
  if ( parmfile.IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
  if (debug>0) mprintf("      Mol2 Title: [%s]\n",buffer);
  //   Scan # atoms and bonds
  // num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
  if ( parmfile.IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
  mol2bonds=0;
  sscanf(buffer,"%i %i",&parmOut.natom, &mol2bonds);
  if (debug>0) {
    mprintf("      Mol2 #atoms: %i\n",parmOut.natom);
    mprintf("      Mol2 #bonds: %i\n",mol2bonds);
  }

  // Allocate memory for atom names, types, and charges.
  parmOut.names = new NAME[ parmOut.natom ];
  parmOut.types = new NAME[ parmOut.natom ];
  parmOut.charge = new double[ parmOut.natom ];
  // Allocate space for coords
  parmOut.parmCoords = new double[ parmOut.natom * 3 ];

  // Get @<TRIPOS>ATOM information
  if (Mol2ScanTo(&parmfile, ATOM)) return 1;
  for (int atom=0; atom < parmOut.natom; atom++) {
    if ( parmfile.IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
    // atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
    //sscanf(buffer,"%*i %s %*f %*f %*f %s %i %s %lf", names[atom], types[atom],
    //       &resnum,resName, charge+atom);
    //mprintf("      %i %s %s %i %s %lf\n",atom,names[atom],types[atom],resnum,resName,charge[atom]);
    Mol2AtomName(buffer,parmOut.names[atom]);
    Mol2AtomType(buffer,parmOut.types[atom]);
    Mol2XYZ(buffer,parmOut.parmCoords + crdidx);
    crdidx += 3;
    Mol2ResNumName(buffer,&resnum,resName);
    parmOut.charge[atom]=Mol2Charge(buffer);
    // Check if residue number has changed - if so record it
    if (resnum != currentResnum) {
      NAME *temprname = new NAME[ parmOut.nres+1 ];
      memcpy(temprname, parmOut.resnames, parmOut.nres * sizeof(NAME));
      delete[] parmOut.resnames;
      parmOut.resnames = temprname;
      strcpy(parmOut.resnames[parmOut.nres], resName);
      int *temprnum = new int[ parmOut.nres+1 ];
      memcpy(temprnum, parmOut.resnums, parmOut.nres * sizeof(int));
      delete[] parmOut.resnums;
      parmOut.resnums = temprnum;
      parmOut.resnums[parmOut.nres]=atom; 
      currentResnum = resnum;
      parmOut.nres++;
    }
  }

  // Get @<TRIPOS>BOND information [optional]
  parmOut.NbondsWithoutH=0;
  parmOut.NbondsWithH=0;
  if (Mol2ScanTo(&parmfile, BOND)==0) {
    for (int bond=0; bond < mol2bonds; bond++) {
      if ( parmfile.IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
      // bond_id origin_atom_id target_atom_id bond_type [status_bits]
      //         resnum         currentResnum
      sscanf(buffer,"%*i %i %i\n",&resnum,&currentResnum);
      // mol2 atom #s start from 1
      parmOut.AddBond(resnum-1, currentResnum-1,0);
    }
  } else {
    mprintf("      Mol2 file does not contain bond information.\n");
  }

  // No box
  parmOut.boxType = NOBOX;

  mprintf("    Mol2 contains %i atoms, %i residues,\n", parmOut.natom,parmOut.nres);
  mprintf("    %i bonds to H, %i other bonds.\n", parmOut.NbondsWithH,parmOut.NbondsWithoutH);

  return 0;
}

