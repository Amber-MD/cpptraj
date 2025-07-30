// Traj_CIF
#include "Traj_CIF.h"
#include "Topology.h"
#include "Frame.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // convertToInteger, convertToDouble

bool Traj_CIF::ID_TrajFormat(CpptrajFile& fileIn) {
  return CIFfile::ID_CIF( fileIn );
}

// Traj_CIF::openTrajin()
// NOTE: Currently CIF files are always read in and stored in memory.
//       No write. Everything handled by setupTrajin and readFrame.
int Traj_CIF::openTrajin() {
  return 0;
}

/** Determine what block has the coordinates we need. */
int Traj_CIF::determineCoordsBlock() {
  blockName_.clear();
  CIFfile::DataBlock const& block1 = file_.GetDataBlock("_atom_site");
  if (!block1.empty()) {
    blockName_ = "_atom_site";
    entryName_ = "_entry";
    xstr_ = "Cartn_x";
    ystr_ = "Cartn_y";
    zstr_ = "Cartn_z";
    nstr_ = "id";
    return 0;
  }
  CIFfile::DataBlock const& block2 = file_.GetDataBlock("_chem_comp_atom");
  if (!block2.empty()) {
    blockName_ = "_chem_comp_atom";
    entryName_ = "_chem_comp";
    xstr_ = "model_Cartn_x";
    ystr_ = "model_Cartn_y";
    zstr_ = "model_Cartn_z";
    nstr_ = "pdbx_ordinal";
    return 0;
  }
  return 1;
}

// Traj_CIF::setupTrajin()
/** Read in entire CIF file. */
int Traj_CIF::setupTrajin(FileName const& fname, Topology* trajParm)
{
  if (file_.Read( fname, debug_ )) return TRAJIN_ERR;
  determineCoordsBlock();
  CIFfile::DataBlock const& block = file_.GetDataBlock(blockName_);
  if (block.empty()) return TRAJIN_ERR;
  // Get coordinate x/y/z columns
  Cartn_x_col_ = block.ColumnIndex(xstr_);
  Cartn_y_col_ = block.ColumnIndex(ystr_);
  Cartn_z_col_ = block.ColumnIndex(zstr_);
  if (Cartn_x_col_ == -1 || Cartn_y_col_ == -1 || Cartn_z_col_ == -1) {
    mprinterr("Error: Could not find %s|%s|%s columns in CIF file.\n",
              xstr_.c_str(), ystr_.c_str(), zstr_.c_str());
    return TRAJIN_ERR;
  }
  // Determine # atoms and # models
  Nmodels_ = 0;
  int model_col = block.ColumnIndex("pdbx_PDB_model_num");
  int id_col    = block.ColumnIndex(nstr_);
  if (id_col == -1) {
    mprinterr("Error: No %s column found in %s block.\n", nstr_.c_str(), blockName_.c_str());
    return TRAJIN_ERR;
  }
  CIFfile::DataBlock::data_it line = block.end();
  --line; // Go to last _atom_site line
  // totalAtoms will == #models * natom
  int totalAtoms = convertToInteger( (*line)[id_col] );
  if (model_col == -1) {
    // No model # column; assume 1 model
    Nmodels_ = 1;
  } else {
    Nmodels_ = convertToInteger( (*line)[model_col] );
  }
  if ( (totalAtoms % Nmodels_) != 0 ) {
    mprintf("Warning: Total number of atoms in CIF (%i) is not divisible by\n"
            "Warning:  number of models (%i). This indicates the number of atoms\n"
            "Warning:  in each model is not the same. Only reading %i atoms of\n"
            "Warning:  the first model.\n",
            totalAtoms, Nmodels_, trajParm->Natom());
    Nmodels_ = 1;
    Natoms_ = trajParm->Natom();
  } else {
    Natoms_ = totalAtoms / Nmodels_;
    if (Natoms_ != trajParm->Natom()) {
      mprinterr("Error: Number of atoms in CIF (%i) does not equal number of atoms\n"
                "Error: in associated topology '%s' (%i)\n", Natoms_,
                trajParm->c_str(), trajParm->Natom());
      return TRAJIN_ERR;
    }
  }
  mprintf("\t%i atoms, %i models.\n", Natoms_, Nmodels_);
  // Get unit cell parameters if present.
  boxInfo_.SetNoBox();
  double cif_box[6];
  int box_stat = file_.cif_Box_verbose( cif_box );
  if (box_stat != -1) {
    if (box_stat || boxInfo_.SetupFromXyzAbg( cif_box )) {
      mprintf("Warning: Box information in CIF appears invalid; disabling box.\n");
      boxInfo_.SetNoBox();
    }
  }
  // Set traj info - No velocity, temperature, time.
  SetCoordInfo( CoordinateInfo( boxInfo_, false, false, false ) );
  // Get title. 
  CIFfile::DataBlock const& entryblock = file_.GetDataBlock(entryName_);
  if (!entryblock.empty())
    SetTitle( entryblock.Data("id") );

  return Nmodels_;
}

// Traj_CIF::readFrame()
int Traj_CIF::readFrame(int set, Frame& frameIn) {
  //if (set >= Nmodels_) return 1;
  // FIXME: Shouldnt have to always search for the block
  CIFfile::DataBlock const& block = file_.GetDataBlock(blockName_);
  CIFfile::DataBlock::data_it line = block.begin() + (set * Natoms_);
  CIFfile::DataBlock::data_it end  = line + Natoms_;
  double *Xptr = frameIn.xAddress(); 
  for (; line != end; ++line) {
    *(Xptr++) = convertToDouble( (*line)[Cartn_x_col_] );
    *(Xptr++) = convertToDouble( (*line)[Cartn_y_col_] );
    *(Xptr++) = convertToDouble( (*line)[Cartn_z_col_] );
  }
  frameIn.SetBox( boxInfo_ );
  return 0;
}

// Traj_CIF::info()
void Traj_CIF::Info() {
  mprintf("is a CIF file");
}
