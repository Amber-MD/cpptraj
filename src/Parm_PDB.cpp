#include "Parm_PDB.h"

int Parm_PDB::ReadParm(std::string const& fname, Topology &TopIn) {
  CpptrajFile infile;
  if (infile.OpenRead(fname)) return 1;
  // Loop over PDB records 
  while ( PDB_GetNextRecord( infile.IOptr() ) ) {
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
  TopIn.SetParmName( pdbtitle, infile.BaseFileStr() );

  infile.CloseFile();
  return 0;
}


bool Parm_PDB::ID_ParmFormat(CpptrajFile& fileIn) {
  // Assumes already set up
  if (fileIn.OpenFile()) return false;
  bool ispdbfile = ID( fileIn.IOptr() );
  fileIn.CloseFile();
  return ispdbfile;
}


