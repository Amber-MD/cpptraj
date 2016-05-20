#include "Parm_PDB.h"
#include "PDBfile.h"
#include "CpptrajStdio.h"
#include "BondSearch.h"
#ifdef TIMER
# include "Timer.h"
#endif

void Parm_PDB::ReadHelp() {
  mprintf("\tpqr     : Read atomic charge/radius from occupancy/B-factor columns (PQR).\n"
          "\treadbox : Read unit cell information from CRYST1 record if present.\n"
          "\tnoconect: Do not read CONECT records if present.\n");
}

int Parm_PDB::processReadArgs(ArgList& argIn) {
  readAsPQR_ = argIn.hasKey("pqr");
  readBox_ = argIn.hasKey("readbox");
  readConect_ = !argIn.hasKey("noconect");
  return 0;
} 

int Parm_PDB::ReadParm(FileName const& fname, Topology &TopIn) {
  PDBfile infile;
  double XYZ[6]; // Hold XYZ/box coords.
  float occupancy, bfactor; // Read in occ/bfac
  BondArray bonds;              // Hold bonds
  std::vector<int> serial;      // Map ATOM/HETATM serial number to actual atom number.
  int atnum;                    // Read in ATOM/HETATM serial number.
  int barray[5];                // Hold CONECT atom and bonds
  char altLoc = ' ';            // For reading in altLoc.
  Frame Coords;
  if (infile.OpenRead(fname)) return 1;
  if (readAsPQR_)   mprintf("\tReading as PQR file.\n");
  if (readBox_)     mprintf("\tUnit cell info will be read from any CRYST1 record.\n");
  if (!readConect_) mprintf("\tNot reading bond info from CONECT records.\n");
# ifdef TIMER
  Timer time_total, time_atom;
  time_total.Start();
# endif
  // Loop over PDB records
  while ( infile.NextRecord() != PDBfile::END_OF_FILE ) {
    if (readBox_ && infile.RecType() == PDBfile::CRYST1) {
      // Box info from CRYST1 record.
      infile.pdb_Box( XYZ );
      TopIn.SetParmBox( XYZ );
    } else if (infile.RecType() == PDBfile::CONECT && readConect_) {
      // BOND - first element will be atom, next few are bonded atoms.
      // To avoid duplicates only add the bond if atom2 > atom1
      int nscan = infile.pdb_Bonds(barray);
      if (nscan > 1) {
        if (nscan > 5) nscan = 5;
        for (int i = 1; i != nscan; i++)
          if (barray[i] > barray[0])
            bonds.push_back( BondType(barray[0], barray[i], -1) );
      }
    } else if (infile.RecType() == PDBfile::ATOM) {
#     ifdef TIMER
      time_atom.Start();
#     endif
      // If this is an ATOM / HETATM keyword, add to topology.
      infile.pdb_XYZ( XYZ );
      Atom pdbAtom = infile.pdb_Atom(altLoc, atnum);
      if (atnum >= (int)serial.size())
        serial.resize( atnum+1, -1 );
      serial[atnum] = TopIn.Natom();
      infile.pdb_OccupancyAndBfactor(occupancy, bfactor);
      if (readAsPQR_) {
        pdbAtom.SetCharge( occupancy );
        pdbAtom.SetGBradius( bfactor );
      } else
        TopIn.AddExtraAtomInfo( AtomExtra(occupancy, bfactor, altLoc) );
      TopIn.AddTopAtom(pdbAtom, infile.pdb_Residue());
      Coords.AddXYZ( XYZ );
#     ifdef TIMER
      time_atom.Stop();
#     endif
    } else if ( infile.RecType() == PDBfile::TER || 
                infile.RecType() == PDBfile::END )
    {
      // Indicate end of molecule for TER/END. Finish if END.
      TopIn.StartNewMol();
      if (infile.RecType() == PDBfile::END) break;
    }
  }
  // Add bonds. The bonds array actually contains ATOM/HETATM serial #s.
  for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
    TopIn.AddBond( serial[bnd->A1()], serial[bnd->A2()] );
  BondSearch( TopIn, Coords, Offset_, debug_ ); 
  // If Topology name not set with TITLE etc, use base filename.
  // TODO: Read in title.
  std::string pdbtitle;
  TopIn.SetParmName( pdbtitle, infile.Filename() );

  infile.CloseFile();
# ifdef TIMER
  time_total.Stop();
  time_atom.WriteTiming(2, "ATOM/HETATM read", time_total.Total());
  time_total.WriteTiming(1, "Total PDB read");
# endif
  return 0;
}

bool Parm_PDB::ID_ParmFormat(CpptrajFile& fileIn) {
  return PDBfile::ID_PDB( fileIn );
}
