#include "Parm_CIF.h"
#include "CIFfile.h"
#include "StringRoutines.h"
#include "CpptrajStdio.h"
#include "BondSearch.h"

/** CONSTRUCTOR */
Parm_CIF::Parm_CIF() :
  read_struct_conn_(false) // TODO enable
{}

// NOTE: MUST correspond to EntryType!
const char* Parm_CIF::Entries_[] = {
  "label_atom_id", "label_comp_id", "Cartn_x", "Cartn_y", "Cartn_z", 
  "label_seq_id", "label_asym_id"
};

/// \return column id of data item with given key; print error message if missing
static inline int checkForCol(CIFfile::DataBlock const& block, std::string const& key)
{
  int colid = block.ColumnIndex( key );
  if (colid == -1)
    mprinterr("Error: CIF block is missing data item '%s'\n", key.c_str());
  return colid;
}

/// Attempt to add a bond between specified atoms to the Topology
static inline void add_cif_bond(Topology& topIn,
                                std::string const& r1name,
                                std::string const& r1num,
                                std::string const& r1chain,
                                std::string const& r1atom,
                                std::string const& r2name,
                                std::string const& r2num,
                                std::string const& r2chain,
                                std::string const& r2atom)
{
  AtomMask m1, m2;
  if (m1.SetMaskString(":" + r1name + "&:;" + r1num + "&::" + r1chain + "&@" + r1atom)) {
    mprinterr("Internal Error: add_cif_bond: Could not set atom 1 mask.\n");
    return;
  }
  if (m2.SetMaskString(":" + r2name + "&:;" + r2num + "&::" + r2chain + "&@" + r2atom)) {
    mprinterr("Internal Error: add_cif_bond: Could not set atom 2 mask.\n");
    return;
  }
  if (topIn.SetupIntegerMask(m1)) {
    mprinterr("Internal Error: add_cif_bond: Could not set up atom 1 mask.\n");
    return;
  }
  if (topIn.SetupIntegerMask(m2)) {
    mprinterr("Internal Error: add_cif_bond: Could not set up atom 2 mask.\n");
    return;
  }
  if (m1.Nselected() > 1) {
    mprinterr("Internal Error: add_cif_bond: atom 1 mask selects more than 1 atom.\n");
    return;
  }
  if (m2.Nselected() > 1) {
    mprinterr("Internal Error: add_cif_bond: atom 2 mask selects more than 1 atom.\n");
    return;
  }
  if (m1.None()) {
    mprintf("Warning: CIF atom 1 '%s' not found.\n", m1.MaskString());
    return;
  }
  if (m2.None()) {
    mprintf("Warning: CIF atom 2 '%s' not found.\n", m2.MaskString());
    return;
  }
}

/** Read data from CIF file. */
int Parm_CIF::ReadParm(FileName const& fname, Topology &TopIn) {
  CIFfile infile;
  if (infile.Read( fname, debug_ )) return 1;

  int read_stat = read_atom_site(infile, fname, TopIn);
  // -1 means the _atom_site data block wasnt found. Try _chem_comp_atom
  if (read_stat == -1)
    read_stat = read_chem_comp_atom(infile, fname, TopIn);
  // If still -1 here then no relevant data blocks found
  if (read_stat == -1) {
    mprinterr("Error: CIF data blocks '_atom_site' or '_chem_comp_atom' not found.\n");
    return 1;
  }

  return read_stat;
}

/// Get essential column, print error message if not found
static inline int get_essential(CIFfile::DataBlock const& block, std::string const& key, FileName const& fname)
{
  int col = block.ColumnIndex( key );
  if (col == -1) {
    mprinterr("Error: In CIF file '%s' could not find entry '%s' in block '%s'\n",
              fname.full(), key.c_str(), block.Header().c_str());
  }
  return col;
}

/** Read _chem_comp_atom data block. */
int Parm_CIF::read_chem_comp_atom(CIFfile& infile, FileName const& fname, Topology& TopIn) const {
  CIFfile::DataBlock const& block = infile.GetDataBlock("_chem_comp_atom");
  if (block.empty()) {
    return -1;
  }
  if (debug_ > 0)
    mprintf("DEBUG: Reading _chem_comp_atom\n");
  // Get essential columns 
  int atom_id_col = get_essential(block, "atom_id", fname);
  if (atom_id_col == -1) return 1;
  int comp_id_col = get_essential(block, "comp_id", fname);
  if (comp_id_col == -1) return 1;
  int model_Cartn_x_col = get_essential(block, "model_Cartn_x", fname);
  if (model_Cartn_x_col == -1) return 1;
  int model_Cartn_y_col = get_essential(block, "model_Cartn_y", fname);
  if (model_Cartn_y_col == -1) return 1;
  int model_Cartn_z_col = get_essential(block, "model_Cartn_z", fname);
  if (model_Cartn_z_col == -1) return 1;
  int type_symbol_col = get_essential(block, "type_symbol", fname);
  if (type_symbol_col == -1) return 1;
  if (debug_ > 0)
    mprintf("DEBUG: cols: %i %i %i %i %i %i\n", atom_id_col, comp_id_col, model_Cartn_x_col, model_Cartn_y_col, model_Cartn_z_col, type_symbol_col);
  // Loop over all atoms
  int current_res = 1; // FIXME only ever one residue?
  double XYZ[3];
  //double occupancy = 1.0;
  //double bfactor = 0.0;
  //char altloc = ' ';
  char icode = ' ';
  std::string chainid = "";
  //int auth_res = -1;
  Frame Coords;
  for (CIFfile::DataBlock::data_it line = block.begin(); line != block.end(); ++line)
  {
    XYZ[0] = convertToDouble( (*line)[ model_Cartn_x_col ] );
    XYZ[1] = convertToDouble( (*line)[ model_Cartn_y_col ] );
    XYZ[2] = convertToDouble( (*line)[ model_Cartn_x_col ] );
    NameType currentResName( (*line)[ comp_id_col ] );
    char ELT[3];
    ELT[0] = ' ';
    ELT[1] = ' ';
    ELT[2] = '\0';
    std::string const& eltstr = (*line)[ type_symbol_col ];
    if (eltstr.size() == 1) {
      ELT[0] = eltstr[0];
    } else if (!eltstr.empty()) {
      ELT[0] = eltstr[0];
      ELT[1] = eltstr[1];
    }
    TopIn.AddTopAtom( Atom((*line)[ atom_id_col ], ELT),
                      Residue(currentResName, current_res, icode, chainid) );
//                              (*line)[ COL[CHAINID] ]) );
    Coords.AddXYZ( XYZ );
  }

  // Determine connectivity
  bool need_bond_search = true;
  //if (read_struct_conn_) {
    CIFfile::DataBlock const& connectBlock = infile.GetDataBlock("_chem_comp_bond");
    if (!connectBlock.empty()) {
      mprintf("\tBlock '_chem_comp_bond' found.\n");
      bool block_is_valid = true;
      int at1_col = get_essential(connectBlock, "atom_id_1", fname);
      int at2_col = get_essential(connectBlock, "atom_id_2", fname);
      block_is_valid = (at1_col != -1 && at2_col != -1);
      if (block_is_valid) {
        for (CIFfile::DataBlock::data_it bond = connectBlock.begin(); bond != connectBlock.end(); ++bond)
        {
          NameType a1name( (*bond)[ at1_col ] );
          NameType a2name( (*bond)[ at2_col ] );
          int idx1 = -1;
          int idx2 = -1;
          for (int atnum = 0; atnum != TopIn.Natom(); atnum++) {
            if ( idx1 == -1 && TopIn[atnum].Name() == a1name ) idx1 = atnum;
            if ( idx2 == -1 && TopIn[atnum].Name() == a2name ) idx2 = atnum;
            if (idx1 > -1 && idx2 > -1) break;
          }
          if (idx1 < 0)
            mprinterr("Error: Atom 1 %s not found for bond to %s\n", *a1name, *a2name);
          if (idx2 < 0)
            mprinterr("Error: Atom 2 %s not found for bond to %s\n", *a2name, *a1name);
          if (idx1 < 0 || idx2 < 0) return 1;
          if (debug_ > 1) 
            mprintf("DEBUG: Adding bond %s %s\n", *a1name, *a2name);
          TopIn.AddBond(idx1, idx2);
        }
        need_bond_search = false;
      }
    } else {
      mprinterr("Error: Required block '_chem_comp_bond' not found.\n");
      return 1;
    }
  //}

  if (need_bond_search) {
    // Search for bonds
    BondSearch bondSearch;
    bondSearch.FindBonds( TopIn, searchType_, Coords, Offset_, debug_ );
  }

  // Get title. 
  CIFfile::DataBlock const& entryblock = infile.GetDataBlock("_chem_comp");
  std::string ciftitle;
  if (!entryblock.empty())
    ciftitle = entryblock.Data("id");
  TopIn.SetParmName( ciftitle, infile.CIFname() );
  // Get unit cell parameters if present.
//  double cif_box[6];
//  int box_stat = infile.cif_Box_verbose( cif_box );
//  if (box_stat != -1) {
//    Box parmBox;
//    parmBox.SetupFromXyzAbg( cif_box );
//    TopIn.SetParmBox( parmBox );
//  }

  return 0;
}

/** Read _atom_site data block. */
int Parm_CIF::read_atom_site(CIFfile& infile, FileName const& fname, Topology& TopIn) const {
  CIFfile::DataBlock::data_it line;

  CIFfile::DataBlock const& block = infile.GetDataBlock("_atom_site");
  if (block.empty()) {
    return -1;
  }
  // Does this CIF contain multiple models?
  int Nmodels = 0;
  int model_col = block.ColumnIndex("pdbx_PDB_model_num");
  if (model_col != -1) {
    line = block.end();
    --line;
    Nmodels = convertToInteger( (*line)[model_col] );
    if (Nmodels > 1)
      mprintf("Warning: CIF '%s' contains %i models. Using first model for topology.\n", 
              fname.full(), Nmodels);
  }
  // Get essential columns
  int COL[NENTRY];
  for (int i = 0; i < (int)NENTRY; i++) {
    COL[i] = block.ColumnIndex(Entries_[i]);
    if (COL[i] == -1) {
      mprinterr("Error: In CIF file '%s' could not find entry '%s' in block '%s'\n",
                fname.full(), Entries_[i], block.Header().c_str());
      return 1;
    }
    if (debug_>0) mprintf("DEBUG: '%s' column = %i\n", Entries_[i], COL[i]);
  }
  // Get optional columns
  int auth_res_col = block.ColumnIndex("auth_seq_id");
  int occ_col = block.ColumnIndex("occupancy");
  int bfac_col = block.ColumnIndex("B_iso_or_equiv");
  int icode_col = block.ColumnIndex("pdbx_PDB_ins_code");
  int altloc_col = block.ColumnIndex("label_alt_id");

  // Loop over all atom sites
  int current_res = 0;
  double XYZ[3];
  double occupancy = 1.0;
  double bfactor = 0.0;
  char altloc = ' ';
  char icode = ' ';
  int auth_res = -1;
  Frame Coords;
  for (line = block.begin(); line != block.end(); ++line) {
    // If more than 1 model check if we are done.
    if (Nmodels > 1) {
      if ( convertToInteger( (*line)[model_col] ) > 1 )
        break;
    }
    if (occ_col != -1) occupancy = convertToDouble( (*line)[ occ_col ] );
    if (bfac_col != -1) bfactor = convertToDouble( (*line)[ bfac_col ] );
    if (altloc_col != -1) altloc = (*line)[ altloc_col ][0];
    // If the 'auth_seq_id' column is present it seems to have residue numbers
    // that mirror PDB residue numbers more closely.
    if (auth_res_col != -1) {
      if (validInteger( (*line)[ auth_res_col ] ))
        auth_res = convertToInteger( (*line)[ auth_res_col ] );
      else
        auth_res = -1;
    }
    // '.' altloc means blank?
    if (altloc == '.') altloc = ' ';
    TopIn.AddAtomAltLoc( altloc );
    TopIn.AddOccupancy( occupancy );
    TopIn.AddBfactor( bfactor );
    if (icode_col != -1) {
      icode = (*line)[ icode_col ][0];
      // '?' icode means blank
      if (icode == '?') icode = ' ';
    }
    XYZ[0] = convertToDouble( (*line)[ COL[X] ] );
    XYZ[1] = convertToDouble( (*line)[ COL[Y] ] );
    XYZ[2] = convertToDouble( (*line)[ COL[Z] ] );
    NameType currentResName( (*line)[ COL[RNAME] ] );
    if ( auth_res != -1 )
      current_res = auth_res;
    else if ( (*line)[ COL[RNUM] ][0] == '.' ) {
      // It seems that in some CIF files, there doesnt have to be a residue
      // number. Check if residue name has changed.
      Topology::res_iterator lastResidue = TopIn.ResEnd() - 1;
      if ( currentResName != lastResidue->Name() )
        current_res = TopIn.Nres() + 1;
    } else
      current_res = convertToInteger( (*line)[ COL[RNUM] ] );
    TopIn.AddTopAtom( Atom((*line)[ COL[ANAME] ], "  "),
                      Residue(currentResName, current_res, icode,
                              (*line)[ COL[CHAINID] ]) );
    Coords.AddXYZ( XYZ );
  }
  // Search for bonds
  BondSearch bondSearch;
  bondSearch.FindBonds( TopIn, searchType_, Coords, Offset_, debug_ );

  if (read_struct_conn_) {
    CIFfile::DataBlock const& connectBlock = infile.GetDataBlock("_struct_conn");
    if (!connectBlock.empty()) {
      mprintf("\tBlock 'struct_conn' found.\n");
      int conn_type_idcol = checkForCol(connectBlock, "conn_type_id");
      int r1_chaincol = checkForCol(connectBlock, "ptnr1_label_asym_id");
      int r1_atomcol  = checkForCol(connectBlock, "ptnr1_label_atom_id");
      int r1_namecol  = checkForCol(connectBlock, "ptnr1_label_comp_id");
      int r1_numcol   = checkForCol(connectBlock, "ptnr1_label_seq_id");
      int r2_chaincol = checkForCol(connectBlock, "ptnr2_label_asym_id");
      int r2_atomcol  = checkForCol(connectBlock, "ptnr2_label_atom_id");
      int r2_namecol  = checkForCol(connectBlock, "ptnr2_label_comp_id");
      int r2_numcol   = checkForCol(connectBlock, "ptnr2_label_seq_id");
      bool block_is_valid = !(conn_type_idcol == -1 ||
                              r1_chaincol == -1 ||
                              r1_atomcol == -1 ||
                              r1_namecol == -1 ||
                              r1_numcol == -1 ||
                              r2_chaincol == -1 ||
                              r2_atomcol == -1 ||
                              r2_namecol == -1 ||
                              r2_numcol == -1
                             );
      if (block_is_valid) {
        // All required columns are present. Loop over struct_conn entries
        for (line = connectBlock.begin(); line != connectBlock.end(); ++line) {
          // Allowed values: covale, disulf, hydrog, metalc
          if ((*line)[conn_type_idcol] == "covale" ||
              (*line)[conn_type_idcol] == "disulf")
          {
            mprintf("%s R1 %s %s Chain ID %s Atom %s -- R2 %s %s Chain ID %s Atom %s\n",
                    (*line)[conn_type_idcol].c_str(),
                    (*line)[r1_namecol].c_str(),
                    (*line)[r1_numcol].c_str(),
                    (*line)[r1_chaincol].c_str(),
                    (*line)[r1_atomcol].c_str(),
                    (*line)[r2_namecol].c_str(),
                    (*line)[r2_numcol].c_str(),
                    (*line)[r2_chaincol].c_str(),
                    (*line)[r2_atomcol].c_str()
                   );
            add_cif_bond(TopIn,
                         (*line)[r1_namecol], (*line)[r1_numcol], (*line)[r1_chaincol], (*line)[r1_atomcol],
                         (*line)[r2_namecol], (*line)[r2_numcol], (*line)[r2_chaincol], (*line)[r2_atomcol]);
          }
        }
      } else {
        mprintf("Warning: Structure connectivity block is missing 1 or more data items, skipping.\n");
      }
    }
  } // END read struct_conn

  // Get title. 
  CIFfile::DataBlock const& entryblock = infile.GetDataBlock("_entry");
  std::string ciftitle;
  if (!entryblock.empty())
    ciftitle = entryblock.Data("id");
  TopIn.SetParmName( ciftitle, infile.CIFname() );
  // Get unit cell parameters if present.
  double cif_box[6];
  int box_stat = infile.cif_Box_verbose( cif_box );
  if (box_stat != -1) {
    Box parmBox;
    parmBox.SetupFromXyzAbg( cif_box );
    TopIn.SetParmBox( parmBox );
  }

  return 0;
}

// Parm_CIF::ID_ParmFormat()
bool Parm_CIF::ID_ParmFormat(CpptrajFile& fileIn) {
  return CIFfile::ID_CIF( fileIn );
}
