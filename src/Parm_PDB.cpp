#include "Parm_PDB.h"

int Parm_PDB::ReadParm(Topology &TopIn) {
  if (OpenFile()) return 1;
  // Indicate parm has coordinates
  TopIn.SetHasCoordinates();
  // Loop over PDB records 
  while ( PDB_GetNextRecord( IO ) ) {
    if (IsPDBatomKeyword()) {
      // If this is an ATOM / HETATM keyword, add to topology
      TopIn.AddAtom(pdb_Atom(), pdb_Residue());
    } else if (IsPDB_TER() || IsPDB_END()) {
      // Indicate end of molecule for TER/END. Finish if END.
      TopIn.StartNewMol();
      if (IsPDB_END()) break;
    }
  }
  // If Topology name not set with TITLE etc, use base filename.
  // TODO: Read in title.
  TopIn.SetParmName( BaseName() );

  CloseFile();
  return 0;
}


bool Parm_PDB::ID_ParmFormat() {
  // Assumes already set up
  if (OpenFile()) return false;
  if (!PDB_GetNextRecord(IO)) return false;
  if (IsPDBkeyword()) return true;
  if (!PDB_GetNextRecord(IO)) return false;
  if (IsPDBkeyword()) return true;
  CloseFile();
  return false;
}


