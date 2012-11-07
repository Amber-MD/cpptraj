#include "Parm_PDB.h"
#include "PDBfile.h"

int Parm_PDB::ReadParm(std::string const& fname, Topology &TopIn) {
  PDBfile infile;
  double XYZ[3];
  if (infile.OpenRead(fname)) return 1;
  // Loop over PDB records 
  while ( infile.NextLine() != 0 ) {
    if (infile.IsPDBatomKeyword()) {
      // If this is an ATOM / HETATM keyword, add to topology
      infile.pdb_XYZ(XYZ);
      TopIn.AddAtom(infile.pdb_Atom(), infile.pdb_Residue(), XYZ);
    } else if (infile.IsPDB_TER() || infile.IsPDB_END()) {
      // Indicate end of molecule for TER/END. Finish if END.
      TopIn.StartNewMol();
      if (infile.IsPDB_END()) break;
    }
  }
  // If Topology name not set with TITLE etc, use base filename.
  // TODO: Read in title.
  std::string pdbtitle;
  TopIn.SetParmName( pdbtitle, infile.BaseFileName() );

  infile.CloseFile();
  return 0;
}

bool Parm_PDB::ID_ParmFormat(CpptrajFile& fileIn) {
  return PDBfile::ID_PDB( fileIn );
}
