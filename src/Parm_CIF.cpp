#include "Parm_CIF.h"
#include "CIFfile.h"
#include "StringRoutines.h"
#include "CpptrajStdio.h"

// NOTE: MUST correspond to EntryType!
const char* Parm_CIF::Entries[] = {
  "label_atom_id", "label_comp_id", "Cartn_x", "Cartn_y", "Cartn_z", 
  "label_seq_id", "label_asym_id"
};

// Parm_CIF::ReadParm()
int Parm_CIF::ReadParm(std::string const& fname, Topology &TopIn) {
  CIFfile infile;
  CIFfile::DataBlock::data_it line;

  if (infile.Read( fname, debug_ )) return 1;
  CIFfile::DataBlock const& block = infile.GetDataBlock("_atom_site");
  if (block.empty()) return 1;
  // Does this CIF contain multiple models?
  int Nmodels = 0;
  int model_col = block.ColumnIndex("pdbx_PDB_model_num");
  if (model_col != -1) {
    line = block.end();
    --line;
    Nmodels = convertToInteger( (*line)[model_col] );
    if (Nmodels > 1)
      mprintf("Warning: CIF '%s' contains %i models. Using first model for topology.\n", 
              fname.c_str(), Nmodels);
  }
  // Get essential columns
  int COL[NENTRY];
  for (int i = 0; i < (int)NENTRY; i++) {
    COL[i] = block.ColumnIndex(Entries[i]);
    if (COL[i] == -1) {
      mprinterr("Error: In CIF file '%s' could not find entry '%s' in block '%s'\n",
                fname.c_str(), Entries[i], block.Header().c_str());
      return 1;
    }
    if (debug_>0) mprintf("DEBUG: '%s' column = %i\n", Entries[i], COL[i]);
  }

  // Loop over all atom sites
  int current_res = 0;
  int last_res = -1;
  double XYZ[3];
  for (line = block.begin(); line != block.end(); ++line) {
    // If more than 1 model check if we are done.
    if (Nmodels > 1) {
      if ( convertToInteger( (*line)[model_col] ) > 1 )
        break;
    }
    XYZ[0] = convertToDouble( (*line)[ COL[X] ] );
    XYZ[1] = convertToDouble( (*line)[ COL[Y] ] );
    XYZ[2] = convertToDouble( (*line)[ COL[Z] ] );
    NameType currentResName( (*line)[ COL[RNAME] ] );
    // It seems that in some CIF files, there doesnt have to be a residue
    // number. Check if residue name has changed.
    if ( (*line)[ COL[RNUM] ][0] == '.' ) {
      Topology::res_iterator lastResidue = TopIn.ResEnd();
      --lastResidue;
      if ( currentResName != (*lastResidue).Name() )
        current_res = TopIn.Nres() + 1;
    } else
      current_res = convertToInteger( (*line)[ COL[RNUM] ] );
    TopIn.AddTopAtom( Atom((*line)[ COL[ANAME] ], (*line)[ COL[CHAINID] ][0], "  "),
                      currentResName, current_res, last_res, XYZ );
  }
  // Get title. 
  CIFfile::DataBlock const& entryblock = infile.GetDataBlock("_entry");
  std::string ciftitle;
  if (!entryblock.empty())
    ciftitle = entryblock.Data("id");
  TopIn.SetParmName( ciftitle, infile.CIFname() );
  
  return 0;
}

// Parm_CIF::ID_ParmFormat()
bool Parm_CIF::ID_ParmFormat(CpptrajFile& fileIn) {
  return CIFfile::ID_CIF( fileIn );
}
