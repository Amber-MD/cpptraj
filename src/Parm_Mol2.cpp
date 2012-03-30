// Parm_Mol2.cpp
#include "Parm_Mol2.h"
#include "CpptrajStdio.h"
#include <cstdio> // sscanf
#include <cstring> // memcpy, strcpy

// Parm_Mol2::ID_ParmFormat() 
bool Parm_Mol2::ID_ParmFormat() {
  // Read the first 10 lines
  if (OpenFile()) return false;
  for (int line = 0; line < 10; line++) {
    if ( IO->Gets(buffer_,BUF_SIZE_) ) return false;
    if ( IsMol2Keyword() ) {
      CloseFile();
      return true;
    }
  }
  CloseFile();
  return false;
}
    

// Parm_Mol2::ReadParm()
/** Read file as a Tripos Mol2 file. */
int Parm_Mol2::ReadParm(Topology &parmOut) {
  if (OpenFile()) return 1;

  mprintf("    Reading Mol2 file %s as topology file.\n",parmOut.c_str());
  // Get @<TRIPOS>MOLECULE information
  if (ScanTo(IO, MOLECULE)) return 1;
  //   Scan title
  if ( IO->Gets(buffer_,BUF_SIZE_) ) return 1;
  if (debug_>0) mprintf("      Mol2 Title: [%s]\n",buffer_);
  // TODO: Set parm title
  //   Scan # atoms and bonds
  // num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
  if ( IO->Gets(buffer_,BUF_SIZE_) ) return 1;
  int natom = 0;
  int mol2bonds = 0;
  sscanf(buffer_,"%i %i",&natom, &mol2bonds);
  if (debug_>0) {
    mprintf("      Mol2 #atoms: %i\n",natom);
    mprintf("      Mol2 #bonds: %i\n",mol2bonds);
  }

  // Allocate memory for atom names, types, and charges.
  //parmOut.names = new NAME[ parmOut.natom ];
  //parmOut.types = new NAME[ parmOut.natom ];
  //parmOut.charge = new double[ parmOut.natom ];
  // Allocate space for coords
  //parmOut.parmCoords = new double[ parmOut.natom * 3 ];

  // Get @<TRIPOS>ATOM information
  if (ScanTo(IO, ATOM)) return 1;
  for (int atom=0; atom < natom; atom++) {
    if ( IO->Gets(buffer_,BUF_SIZE_) ) return 1;
    parmOut.AddAtom( Mol2Atom(), Mol2Residue() );
  }

  // Get @<TRIPOS>BOND information [optional]
  int at1 = 0;
  int at2 = 0;
  if (ScanTo(IO, BOND)==0) {
    for (int bond=0; bond < mol2bonds; bond++) {
      if ( IO->Gets(buffer_,BUF_SIZE_) ) return 1;
      // bond_id origin_atom_id target_atom_id bond_type [status_bits]
      //         resnum         currentResnum
      sscanf(buffer_,"%*i %i %i\n", &at1, &at2);
      // mol2 atom #s start from 1
      parmOut.AddBond(at1-1, at2-1);
    }
  } else {
    mprintf("      Mol2 file does not contain bond information.\n");
  }

  // No box
  parmOut.SetNoBox();

  mprintf("    Mol2 contains %i atoms, %i residues,\n", parmOut.Natom(),parmOut.Nres());
  //mprintf("    %i bonds to H, %i other bonds.\n", parmOut.NbondsWithH,parmOut.NbondsWithoutH);

  CloseFile();

  return 0;
}

