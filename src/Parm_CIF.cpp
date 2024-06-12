#include "Parm_CIF.h"
#include "CIFfile.h"
#include "StringRoutines.h"
#include "CpptrajStdio.h"
#include "BondSearch.h"

// NOTE: MUST correspond to EntryType!
const char* Parm_CIF::Entries[] = {
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

// Parm_CIF::ReadParm()
int Parm_CIF::ReadParm(FileName const& fname, Topology &TopIn) {
  CIFfile infile;
  CIFfile::DataBlock::data_it line;

  if (infile.Read( fname, debug_ )) return 1;
  CIFfile::DataBlock const& block = infile.GetDataBlock("_atom_site");
  if (block.empty()) {
    mprinterr("Error: CIF data block '_atom_site' not found.\n");
    return 1;
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
    COL[i] = block.ColumnIndex(Entries[i]);
    if (COL[i] == -1) {
      mprinterr("Error: In CIF file '%s' could not find entry '%s' in block '%s'\n",
                fname.full(), Entries[i], block.Header().c_str());
      return 1;
    }
    if (debug_>0) mprintf("DEBUG: '%s' column = %i\n", Entries[i], COL[i]);
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

  CIFfile::DataBlock const& connectBlock = infile.GetDataBlock("_struct_conn");
  if (!connectBlock.empty()) {
    mprintf("\tBlock 'struct_conn' found.\n");
    int conn_type_idcol = checkForCol(connectBlock, "conn_type_id");
    int r1_chaincol = checkForCol(connectBlock, "ptnr1_label_asym_id");
    int r1_atomcol  = checkForCol(connectBlock, "ptnr1_label_atom_id");
    int r1_namecol  = checkForCol(connectBlock, "ptnr1_label_comp_id");
    int r1_numcol   = checkForCol(connectBlock, "ptnr1_label_seq_id");
    bool block_is_valid = !(conn_type_idcol == -1 ||
                            r1_chaincol == -1 ||
                            r1_atomcol == -1 ||
                            r1_namecol == -1 ||
                            r1_numcol == -1
                           );
    if (block_is_valid) {
      for (line = connectBlock.begin(); line != connectBlock.end(); ++line) {
        mprintf("R1 %s %s Chain ID %s Atom %s\n",
                (*line)[r1_namecol].c_str(),
                (*line)[r1_numcol].c_str(),
                (*line)[r1_chaincol].c_str(),
                (*line)[r1_atomcol].c_str());
      }
    } else {
      mprintf("Warning: Structure connectivity block is missing 1 or more data items, skipping.\n");
    }
  }
  // Get title. 
  CIFfile::DataBlock const& entryblock = infile.GetDataBlock("_entry");
  std::string ciftitle;
  if (!entryblock.empty())
    ciftitle = entryblock.Data("id");
  TopIn.SetParmName( ciftitle, infile.CIFname() );
  // Get unit cell parameters if present.
  CIFfile::DataBlock const& cellblock = infile.GetDataBlock("_cell");
  if (!cellblock.empty()) {
    double cif_box[6];
    cif_box[0] = convertToDouble( cellblock.Data("length_a") );
    cif_box[1] = convertToDouble( cellblock.Data("length_b") );
    cif_box[2] = convertToDouble( cellblock.Data("length_c") );
    cif_box[3] = convertToDouble( cellblock.Data("angle_alpha") );
    cif_box[4] = convertToDouble( cellblock.Data("angle_beta" ) );
    cif_box[5] = convertToDouble( cellblock.Data("angle_gamma") );
    mprintf("\tRead cell info from CIF: a=%g b=%g c=%g alpha=%g beta=%g gamma=%g\n",
              cif_box[0], cif_box[1], cif_box[2], cif_box[3], cif_box[4], cif_box[5]);
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
