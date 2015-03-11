#include "Parm_PDB.h"
#include "PDBfile.h"
#include "CpptrajStdio.h"

void Parm_PDB::ReadHelp() {
  mprintf("\tpqr:     Read atomic charge/radius from occupancy/B-factor columns (PQR).\n"
          "\treadbox: Read unit cell information from CRYST1 record if present.\n");
}

int Parm_PDB::processReadArgs(ArgList& argIn) {
  readAsPQR_ = argIn.hasKey("pqr");
  readBox_ = argIn.hasKey("readbox");
  return 0;
} 

int Parm_PDB::ReadParm(std::string const& fname, Topology &TopIn) {
  PDBfile infile;
  double XYZ[6]; // Hold XYZ/box coords.
  float occupancy, bfactor; // Read in occ/bfac
  std::vector<AtomExtra> extra; // Hold occ/bfac if not PQR
  std::vector<NameType> Icodes; // Hold residue icodes
  char icode[2];                // For reading in icode.
  char altLoc = ' ';            // For reading in altLoc.
  icode[1] = '\0';
  if (infile.OpenRead(fname)) return 1;
  // Loop over PDB records
  while ( infile.NextRecord() != PDBfile::END_OF_FILE ) {
    if (readBox_ && infile.RecType() == PDBfile::CRYST1) {
      // Box info from CRYST1 record.
      infile.pdb_Box( XYZ );
      TopIn.SetParmBox( XYZ );
    } else if (infile.RecType() == PDBfile::ATOM) {
      // If this is an ATOM / HETATM keyword, add to topology.
      infile.pdb_XYZ( XYZ );
      Atom pdbAtom = infile.pdb_Atom(altLoc);
      infile.pdb_OccupanyAndBfactor(occupancy, bfactor);
      if (readAsPQR_) {
        pdbAtom.SetCharge( occupancy );
        pdbAtom.SetGBradius( bfactor );
      } else
        extra.push_back( AtomExtra(occupancy, bfactor, altLoc) );
      TopIn.AddTopAtom(pdbAtom, infile.pdb_ResNum(icode[0]), infile.pdb_ResName(), XYZ);
      Icodes.push_back( NameType(icode) );
    } else if ( infile.RecType() == PDBfile::TER || 
                infile.RecType() == PDBfile::END )
    {
      // Indicate end of molecule for TER/END. Finish if END.
      TopIn.StartNewMol();
      if (infile.RecType() == PDBfile::END) break;
    }
  }
  if (TopIn.SetExtraAtomInfo(0, extra, Icodes)) return 1;
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
