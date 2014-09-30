#include "Parm_PDB.h"
#include "PDBfile.h"
#include "CpptrajStdio.h"

void Parm_PDB::ReadHelp() {
  mprintf("\t[pqr] [readbox]\n");
}

int Parm_PDB::processReadArgs(ArgList& argIn) {
  readAsPQR_ = argIn.hasKey("pqr");
  readBox_ = argIn.hasKey("readbox");
  return 0;
} 

int Parm_PDB::ReadParm(std::string const& fname, Topology &TopIn) {
  PDBfile infile;
  double XYZ[6];
  int current_res = 0;
  if (infile.OpenRead(fname)) return 1;
  // Loop over PDB records 
  while ( infile.NextLine() != 0 ) {
    if (readBox_ && infile.IsBoxKeyword()) {
      // CRYST1 keyword: RECORD A B C ALPHA BETA GAMMA SGROUP Z
      infile.pdb_Box( XYZ );
      mprintf("\tRead CRYST1 info from PDB: a=%g b=%g c=%g alpha=%g beta=%g gamma=%g\n",
              XYZ[0], XYZ[1], XYZ[2], XYZ[3], XYZ[4], XYZ[5]);
      // Warn if the box looks strange.
      if (XYZ[0] == 1.0 && XYZ[0] == XYZ[1] && XYZ[0] == XYZ[2])
        mprintf("Warning: PDB cell lengths are all 1.0 Ang.;"
                " this usually indicates an invalid box.\n");
      TopIn.SetBox( Box(XYZ) );
    } else if (infile.IsPDBatomKeyword()) {
      // If this is an ATOM / HETATM keyword, add to topology
      infile.pdb_XYZ(XYZ);
      NameType pdbresname = infile.pdb_Residue( current_res );
      TopIn.AddTopAtom(infile.pdb_Atom(readAsPQR_), current_res, pdbresname, XYZ);
    } else if (infile.IsPDB_TER() || infile.IsPDB_END()) {
      // Indicate end of molecule for TER/END. Finish if END.
      TopIn.StartNewMol();
      if (infile.IsPDB_END()) break;
    }
  }
  // If Topology name not set with TITLE etc, use base filename.
  // TODO: Read in title.
  std::string pdbtitle;
  TopIn.SetParmName( pdbtitle, infile.Filename() );

  infile.CloseFile();
  return 0;
}

bool Parm_PDB::ID_ParmFormat(CpptrajFile& fileIn) {
  return PDBfile::ID_PDB( fileIn );
}
