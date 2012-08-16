// Parm_Mol2.cpp
#include "Parm_Mol2.h"
#include "CpptrajStdio.h"

// Parm_Mol2::ID_ParmFormat() 
bool Parm_Mol2::ID_ParmFormat() {
  if (OpenFile()) return false;
  //mprintf("DEBUG: Checking Mol2 parm format.\n");
  bool ismol2file = ID( IO );
  CloseFile();
  return ismol2file;
}
    
// Parm_Mol2::ReadParm()
/** Read file as a Tripos Mol2 file. */
int Parm_Mol2::ReadParm(Topology &parmOut) {
  if (OpenFile()) return 1;
  mprintf("    Reading Mol2 file %s as topology file.\n",BaseFileStr());
  // Get @<TRIPOS>MOLECULE information
  if (ReadMolecule(IO)) return 1;
  parmOut.SetParmName( Mol2Title(), BaseFileStr() );

  // Allocate memory for atom names, types, and charges.
  //parmOut.names = new NAME[ parmOut.natom ];
  //parmOut.types = new NAME[ parmOut.natom ];
  //parmOut.charge = new double[ parmOut.natom ];
  // Allocate space for coords
  //parmOut.parmCoords = new double[ parmOut.natom * 3 ];

  // Get @<TRIPOS>ATOM information
  if (ScanTo(IO, ATOM)) return 1;
  for (int atom=0; atom < Mol2Natoms(); atom++) {
    if ( GetLine( IO ) ) return 1;
    parmOut.AddAtom( Mol2Atom(), Mol2Residue(), XYZ() );
  }

  // Get @<TRIPOS>BOND information [optional]
  int at1 = 0;
  int at2 = 0;
  if (ScanTo(IO, BOND)==0) {
    for (int bond=0; bond < Mol2Nbonds(); bond++) {
      if ( GetLine( IO ) ) return 1;
      // bond_id origin_atom_id target_atom_id bond_type [status_bits]
      //         resnum         currentResnum
      Mol2Bond(at1, at2);
      // mol2 atom #s start from 1
      parmOut.AddBond(at1-1, at2-1);
    }
  } else {
    mprintf("      Mol2 file does not contain bond information.\n");
  }

  // No box
  parmOut.ParmBox().SetNoBox();

  mprintf("    Mol2 contains %i atoms, %i residues,\n", parmOut.Natom(),parmOut.Nres());
  //mprintf("    %i bonds to H, %i other bonds.\n", parmOut.NbondsWithH,parmOut.NbondsWithoutH);

  CloseFile();

  return 0;
}

