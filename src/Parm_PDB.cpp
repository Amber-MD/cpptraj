#include "Parm_PDB.h"
#include "PDBfile.h"
#include "CpptrajStdio.h"
#include "BondSearch.h"
#ifdef TIMER
# include "Timer.h"
#endif

void Parm_PDB::ReadHelp() {
  mprintf("\tpqr      : Read atomic charge/radius from occupancy/B-factor columns (PQR).\n"
          "\treadbox  : Read unit cell information from CRYST1 record if present.\n"
          "\tconect   : Read CONECT records if present (default).\n"
          "\tnoconect : Do not read CONECT records if present.\n");
}

int Parm_PDB::processReadArgs(ArgList& argIn) {
  readAsPQR_ = argIn.hasKey("pqr");
  readBox_ = argIn.hasKey("readbox");
  if (argIn.hasKey("conect"))
    ConectMode_ = READ;
  else if (argIn.hasKey("noconect"))
    ConectMode_ = SKIP;
  return 0;
}

int Parm_PDB::ReadParm(FileName const& fname, Topology &TopIn) {
  typedef std::vector<PDBfile::Link> Larray;
  PDBfile infile;
  double XYZ[6]; // Hold XYZ/box coords.
  float occupancy, bfactor; // Read in occ/bfac
  BondArray bonds;              // Hold bonds
  Larray links;                 // Hold LINK bonds
  std::vector<int> serial;      // Map ATOM/HETATM serial number to actual atom number.
  int atnum;                    // Read in ATOM/HETATM serial number.
  int barray[5];                // Hold CONECT atom and bonds
  char altLoc = ' ';            // For reading in altLoc.
  Frame Coords;
  // Determine if CONECT records should be read.
  bool readConect;
  if (ConectMode_ == SKIP)
    readConect = false;
  else
    readConect = true;
  if (infile.OpenRead(fname)) return 1;
  if (readAsPQR_)  mprintf("\tReading as PQR file.\n");
  if (readBox_)    mprintf("\tUnit cell info will be read from any CRYST1 record.\n");
  if (!readConect) mprintf("\tNot reading bond info from CONECT records.\n");
# ifdef TIMER
  Timer time_total, time_atom;
  time_total.Start();
# endif
  bool missingResidues = false;
  // Loop over PDB records
  while ( infile.NextRecord() != PDBfile::END_OF_FILE ) {
    if (readBox_ && infile.RecType() == PDBfile::CRYST1) {
      // Box info from CRYST1 record.
      infile.pdb_Box( XYZ );
      TopIn.SetParmBox( XYZ );
    } else if (infile.RecType() == PDBfile::CONECT && readConect) {
      // BOND - first element will be atom, next few are bonded atoms.
      // To avoid duplicates only add the bond if atom2 > atom1
      int nscan = infile.pdb_Bonds(barray);
      if (nscan > 1) {
        if (nscan > 5) nscan = 5;
        for (int i = 1; i != nscan; i++)
          if (barray[i] > barray[0])
            bonds.push_back( BondType(barray[0], barray[i], -1) );
      }
    } else if (infile.RecType() == PDBfile::LINK && readConect) {
      // LINK
      links.push_back( infile.pdb_Link() );
      if (debug_ > 0) {
        PDBfile::Link const& lr = links.back();
        mprintf("DEBUG: Link record: %s %s %i to %s %s %i\n",
                lr.aname1(), lr.rname1(), lr.Rnum1(),
                lr.aname2(), lr.rname2(), lr.Rnum2());
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
      if (readAsPQR_) {
        infile.pdb_ChargeAndRadius(occupancy, bfactor);
        pdbAtom.SetCharge( occupancy );
        pdbAtom.SetGBradius( bfactor );
      } else {
        infile.pdb_OccupancyAndBfactor(occupancy, bfactor);
        TopIn.AddExtraAtomInfo( AtomExtra(occupancy, bfactor, altLoc) );
      }
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
    } else if ( !missingResidues && infile.RecType() == PDBfile::MISSING_RES ) {
      missingResidues = true;
      mprintf("Warning: PDB file has MISSING RESIDUES section.\n");
      if (ConectMode_ == UNSPECIFIED) {
        mprintf("Warning: Not reading any connectivity (CONECT etc).\n"
                "Warning: Use the 'conect' option to force reading of these records.\n");
        readConect = false;
        links.clear();
        bonds.clear();
      } else if (ConectMode_ == READ)
        mprintf("Warning: If molecule determination fails try specifying 'noconect' instead.\n");
    }
  }
  // Sanity check
  if (TopIn.Natom() < 1) {
    mprinterr("Error: No atoms present in PDB.\n");
    return 1;
  }
  // Add bonds. The bonds array actually contains ATOM/HETATM serial #s.
  if (!bonds.empty()) {
    if (serial.empty()) // Should never get here
      mprintf("Warning: CONECT info present but no ATOM/HETATMs.\n");
    else
      for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
        TopIn.AddBond( serial[bnd->A1()], serial[bnd->A2()] );
  }
  // Add LINK bonds. Need to search for original residue numbers here.
  if (!links.empty()) {
    for (Larray::const_iterator link = links.begin(); link != links.end(); ++link) {
      Topology::res_iterator r1 = TopIn.ResEnd();
      Topology::res_iterator r2 = TopIn.ResEnd();
      for (Topology::res_iterator res = TopIn.ResStart(); res != TopIn.ResEnd(); ++res) {
        if (r1 == TopIn.ResEnd()) {
          if (link->Rnum1() == res->OriginalResNum() &&
              link->Chain1() == res->ChainID() &&
              link->Icode1() == res->Icode())
          {
            r1 = res;
            if (r2 != TopIn.ResEnd()) break;
          }
        }
        if (r2 == TopIn.ResEnd()) {
          if (link->Rnum2() == res->OriginalResNum() &&
              link->Chain2() == res->ChainID() &&
              link->Icode2() == res->Icode())
          {
            r2 = res;
            if (r1 != TopIn.ResEnd()) break;
          }
        }
      } // END loop over topology residues
      // SANITY CHECK
      if (r1 == TopIn.ResEnd()) {
        mprintf("Warning: Could not find 1st residue %i %s for LINK record.\n", link->Rnum1(), link->rname1());
      } else if (r2 == TopIn.ResEnd()) {
        mprintf("Warning: Could not find 2nd residue %i %s for LINK record.\n", link->Rnum2(), link->rname2());
      } else {
        int idx1 = TopIn.FindAtomInResidue(r1 - TopIn.ResStart(), NameType(link->aname1()));
        if (idx1 < 0) {
          mprintf("Warning: Could not find 1st atom %s in residue %i %s for LINK record.\n", link->aname1(), link->Rnum1(), link->rname1());
        } else {
          int idx2 = TopIn.FindAtomInResidue(r2 - TopIn.ResStart(), NameType(link->aname2()));
          if (idx2 < 0) {
            mprintf("Warning: Could not find 2nd atom %s in residue %i %s for LINK record.\n", link->aname2(), link->Rnum2(), link->rname2());
          } else {
            if (debug_ > 0)
              mprintf("DEBUG: Adding bond %s to %s\n",
                      TopIn.TruncResAtomNameNum(idx1).c_str(),
                      TopIn.TruncResAtomNameNum(idx2).c_str());
            TopIn.AddBond(idx1, idx2);
          }
        }
      }
    } // END loop over Link records
  } // END PDB has link records.
  // Fill in bonds
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
