#include "Parm_PDB.h"

int Parm_PDB::ReadParm(Topology &TopIn) {
  if (OpenFile()) return 1;
  // Loop over PDB records 
  while ( PDB_GetNextRecord( IO ) ) {
    if (IsPDBatomKeyword()) {
      // If this is an ATOM / HETATM keyword, add to topology
      TopIn.AddAtom(pdb_Atom(), pdb_Residue(), XYZ());
    } else if (IsPDB_TER() || IsPDB_END()) {
      // Indicate end of molecule for TER/END. Finish if END.
      TopIn.StartNewMol();
      if (IsPDB_END()) break;
    }
  }
  // If Topology name not set with TITLE etc, use base filename.
  // TODO: Read in title.
  std::string pdbtitle;
  TopIn.SetParmName( pdbtitle, BaseFileStr() );

  CloseFile();
  return 0;
}


bool Parm_PDB::ID_ParmFormat() {
  // Assumes already set up
  if (OpenFile()) return false;
  bool ispdbfile = ID( IO );
  CloseFile();
  return ispdbfile;
}


