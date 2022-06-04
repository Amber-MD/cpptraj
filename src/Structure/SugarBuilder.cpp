#include "SugarBuilder.h"
#include "../ArgList.h"
#include "../CpptrajFile.h"
#include "../CpptrajStdio.h"

using namespace Cpptraj::Structure;

/** Load reduced internal PDB to Glycam map. */
void SugarBuilder::SetGlycamPdbResMap() {
  pdb_to_glycam_.insert( PairType("64K",
    SugarToken("alpha-D-arabinopyranose", "A", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("AHR",
    SugarToken("alpha-L-arabinofuranose", "A", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("ARA",
    SugarToken("alpha-L-arabinopyranose", "A", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("ARB",
    SugarToken("beta-L-arabinopyranose", "A", SugarToken::BETA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("AXR",
    SugarToken("methyl alpha-D-arabinofuranoside", "A", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("BXX",
    SugarToken("beta-D-arabinofuranose", "A", SugarToken::BETA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("BXY",
    SugarToken("alpha-D-arabinofuranose", "A", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("FUB",
    SugarToken("beta-L-arabinofuranose", "A", SugarToken::BETA, SugarToken::IS_L, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("SEJ",
    SugarToken("beta-D-arabinopyranose", "A", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("LDY",
    SugarToken("alpha-D-lyxopyranose", "D", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z4W",
    SugarToken("beta-D-lyxopyranose", "D", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("0MK",
    SugarToken("beta-L-ribopyranose", "R", SugarToken::BETA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("32O",
    SugarToken("beta-L-ribofuranose", "R", SugarToken::BETA, SugarToken::IS_L, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("BDR",
    SugarToken("beta-D-ribofuranose", "R", SugarToken::BETA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("RIB",
    SugarToken("alpha-D-ribofuranose", "R", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("RIP",
    SugarToken("beta-D-ribopyranose", "R", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("YYM",
    SugarToken("alpha-D-ribopyranose", "R", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z6J",
    SugarToken("alpha-L-ribofuranose", "R", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("HSY",
    SugarToken("alpha-L-xylopyranose", "X", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("HSZ",
    SugarToken("beta-D-xylopyranose", "X", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("LXC",
    SugarToken("beta-L-xylopyranose", "X", SugarToken::BETA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("XYP",
    SugarToken("beta-D-xylopyranose", "X", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("XYS",
    SugarToken("alpha-D-xylopyranose", "X", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("XYZ",
    SugarToken("beta-D-xylofuranose", "X", SugarToken::BETA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("AFD",
    SugarToken("alpha-D-allopyranose", "N", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("ALL",
    SugarToken("beta-D-allopyranose", "N", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("VDS",
    SugarToken("beta-D-allofuranose", "N", SugarToken::BETA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("VDV",
    SugarToken("alpha-D-allofuranose", "N", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("WOO",
    SugarToken("beta-L-allopyranose", "N", SugarToken::BETA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z2D",
    SugarToken("alpha-L-allopyranose", "N", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("3MK",
    SugarToken("beta-L-altropyranose", "E", SugarToken::BETA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("SHD",
    SugarToken("alpha-D-altropyranose", "E", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z6H",
    SugarToken("alpha-L-altropyranose", "E", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GAL",
    SugarToken("beta-D-galactopyranose", "L", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GIV",
    SugarToken("beta-L-galactopyranose", "L", SugarToken::BETA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GLA",
    SugarToken("alpha-D-galactopyranose", "L", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GXL",
    SugarToken("alpha-L-galactopyranose", "L", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GZL",
    SugarToken("beta-D-galactofuranose", "L", SugarToken::BETA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("BGC",
    SugarToken("beta-D-glucopyranose", "G", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GLC",
    SugarToken("alpha-D-glucopyranose", "G", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GU4",
    SugarToken("2,3,4,6-tetra-O-sulfonato-alpha-D-glucopyranose", "G", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("MGL",
    SugarToken("methyl beta-D-glucopyranoside", "G", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z8T",
    SugarToken("beta-L-glucopyranose", "G", SugarToken::BETA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("4GL",
    SugarToken("alpha-D-gulopyranose", "K", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GL0",
    SugarToken("beta-D-gulopyranose", "K", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GUP",
    SugarToken("alpha-L-gulopyranose", "K", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z8H",
    SugarToken("beta-L-gulopyranose", "K", SugarToken::BETA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("BMA",
    SugarToken("beta-D-mannopyranose", "M", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("MAN",
    SugarToken("alpha-D-mannopyranose", "M", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("A5C",
    SugarToken("alpha-L-talofuranose", "T", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("SDY",
    SugarToken("beta-D-talopyranose", "T", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("ZEE",
    SugarToken("beta-L-talopyranose", "T", SugarToken::BETA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("BDF",
    SugarToken("beta-D-fructopyranose", "C", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("FRU",
    SugarToken("beta-D-fructofuranose", "C", SugarToken::BETA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("LFR",
    SugarToken("beta-L-fructofuranose", "C", SugarToken::BETA, SugarToken::IS_L, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("YYJ",
    SugarToken("1,3,4,6-tetra-O-sulfo-beta-D-fructofuranose", "C", SugarToken::BETA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("Z9N",
    SugarToken("alpha-D-fructofuranose", "C", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("PSV",
    SugarToken("alpha-D-psicofuranose", "P", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("SF6",
    SugarToken("alpha-L-psicofuranose", "P", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("SF9",
    SugarToken("beta-L-psicofuranose", "P", SugarToken::BETA, SugarToken::IS_L, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("TTV",
    SugarToken("beta-D-psicofuranose", "P", SugarToken::BETA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("SOE",
    SugarToken("alpha-L-sorbopyranose", "B", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("UEA",
    SugarToken("beta-D-sorbofuranose", "B", SugarToken::BETA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("T6T",
    SugarToken("alpha-D-tagatopyranose", "J", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("FCA",
    SugarToken("alpha-D-fucopyranose", "F", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("FCB",
    SugarToken("beta-D-fucopyranose", "F", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("FUC",
    SugarToken("alpha-L-fucopyranose", "F", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("FUL",
    SugarToken("beta-L-fucopyranose", "F", SugarToken::BETA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GYE",
    SugarToken("beta-D-fucofuranose", "F", SugarToken::BETA, SugarToken::IS_D, SugarToken::FURANOSE)) );
  pdb_to_glycam_.insert( PairType("MXY",
    SugarToken("2-O-methyl-beta-L-fucopyranose", "F", SugarToken::BETA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("MXZ",
    SugarToken("2-O-methyl-alpha-L-fucopyranose", "F", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("G6D",
    SugarToken("alpha-D-quinovopyranose", "Q", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("YYK",
    SugarToken("beta-D-quinovopyranose", "Q", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("RAM",
    SugarToken("alpha-L-rhamnopyranose", "H", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("RM4",
    SugarToken("beta-L-rhamnopyranose", "H", SugarToken::BETA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("XXR",
    SugarToken("alpha-D-rhamnopyranose", "H", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("ADA",
    SugarToken("alpha-D-galactopyranuronic acid", "O", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GTR",
    SugarToken("beta-D-galactopyranuronic acid", "O", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("BDP",
    SugarToken("beta-D-glucopyranuronic acid", "Z", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GCU",
    SugarToken("alpha-D-glucopyranuronic acid", "Z", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GCV",
    SugarToken("4-O-methyl-alpha-D-glucopyranuronic acid", "Z", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("IDR",
    SugarToken("alpha-L-idopyranuronic acid", "U", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("IDS",
    SugarToken("2-O-sulfo-alpha-L-idopyranuronic acid", "U", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("A2G",
    SugarToken("2-acetamido-2-deoxy-alpha-D-galactopyranose", "V", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("ASG",
    SugarToken("2-acetamido-2-deoxy-4-O-sulfo-beta-D-galactopyranose", "V", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("NG6",
    SugarToken("2-acetamido-2-deoxy-6-O-sulfo-beta-D-galactopyranose", "V", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("NGA",
    SugarToken("2-acetamido-2-deoxy-beta-D-galactopyranose", "V", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("YYQ",
    SugarToken("2-acetamido-2-deoxy-alpha-L-galactopyranose", "V", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("NAG",
    SugarToken("2-acetamido-2-deoxy-beta-D-glucopyranose", "Y", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("NDG",
    SugarToken("2-acetamido-2-deoxy-alpha-D-glucopyranose", "Y", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("NGZ",
    SugarToken("2-acetamido-2-deoxy-alpha-L-glucopyranose", "Y", SugarToken::ALPHA, SugarToken::IS_L, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("BM3",
    SugarToken("2-acetamido-2-deoxy-alpha-D-mannopyranose", "W", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("BM7",
    SugarToken("2-acetamido-2-deoxy-beta-D-mannopyranose", "W", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("SIA",
    SugarToken("N-acetyl-alpha-neuraminic acid", "S", SugarToken::ALPHA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("SLB",
    SugarToken("N-acetyl-beta-neuraminic acid", "S", SugarToken::BETA, SugarToken::IS_D, SugarToken::PYRANOSE)) );
  // PDB to glycam atom name maps
  // 0 - V,W,Y C7,C2N  O7,O2N  C8,CME
  pdb_glycam_name_maps_.push_back(NameMapType());
  pdb_glycam_name_maps_A_.push_back(NameMapType());
  pdb_glycam_name_maps_B_.push_back(NameMapType());
  pdb_glycam_name_maps_.back().insert(NamePairType("C7","C2N"));
  pdb_glycam_name_maps_.back().insert(NamePairType("O7","O2N"));
  pdb_glycam_name_maps_.back().insert(NamePairType("C8","CME"));
  glycam_res_idx_map_.insert( ResIdxPairType( "V", 0 ) );
  glycam_res_idx_map_.insert( ResIdxPairType( "W", 0 ) );
  glycam_res_idx_map_.insert( ResIdxPairType( "Y", 0 ) );
  // 1 - S     C10,C5N O10,O5N C11,CME
  pdb_glycam_name_maps_.push_back(NameMapType());
  pdb_glycam_name_maps_A_.push_back(NameMapType());
  pdb_glycam_name_maps_B_.push_back(NameMapType());
  pdb_glycam_name_maps_.back().insert(NamePairType("C10","C5N"));
  pdb_glycam_name_maps_.back().insert(NamePairType("O10","O5N"));
  pdb_glycam_name_maps_.back().insert(NamePairType("C11","CME"));
  glycam_res_idx_map_.insert( ResIdxPairType( "S", 1 ) );
  // 2 - H     C6,C6M,B
  pdb_glycam_name_maps_.push_back(NameMapType());
  pdb_glycam_name_maps_A_.push_back(NameMapType());
  pdb_glycam_name_maps_B_.push_back(NameMapType());
  pdb_glycam_name_maps_B_.back().insert(NamePairType("C6","C6M"));
  glycam_res_idx_map_.insert( ResIdxPairType( "H", 2 ) );
  // PDB to glycam linkage residue name maps
  pdb_glycam_linkageRes_map_.insert( NamePairType("SER", "OLS") );
  pdb_glycam_linkageRes_map_.insert( NamePairType("THR", "OLT") );
  pdb_glycam_linkageRes_map_.insert( NamePairType("HYP", "OLP") );
  pdb_glycam_linkageRes_map_.insert( NamePairType("ASN", "NLN") );
}

/** Load PDB to Glycam residue map from file. */
int SugarBuilder::LoadGlycamPdbResMap(std::string const& fnameIn)
{
  std::string fname = fnameIn;
  if (fnameIn.empty()) {
    // Check CPPTRAJHOME
    const char* env = getenv("CPPTRAJHOME");
    if (env != 0) {
      fname.assign(env);
      fname += "/dat/Carbohydrate_PDB_Glycam_Names.txt";
      mprintf("Info: Parameter file path from CPPTRAJHOME variable: '%s'\n", fname.c_str());
    } else {
      // Check AMBERHOME
      env = getenv("AMBERHOME");
      if (env != 0) {
        fname.assign(env);
        fname += "/AmberTools/src/cpptraj/dat/Carbohydrate_PDB_Glycam_Names.txt";
        mprintf("Info: Parameter file path from AMBERHOME variable: '%s'\n", fname.c_str());
      }
    }
  }
  if (fname.empty()) {
    mprintf("Warning: No PDB->Glycam file specified and/or CPPTRAJHOME not set.\n"
            "Warning: Using only basic PDB residue name recognition.\n");
    SetGlycamPdbResMap();
    return 0;
  }
  mprintf("\tReading PDB residue name -> Glycam name map from '%s'\n", fname.c_str());

  CpptrajFile infile;
  if (infile.OpenRead(fname)) {
    mprinterr("Error: Could not open Glycam residue map file.\n");
    return 1;
  }
  const char* ptr = 0;
  // Describe which section of the file we are in
  enum SectionType { PDB_RESMAP_SECTION = 0, PDB_ATOMMAP_SECTION, PDB_LINKAGE_RES_SECTION };
  SectionType section = PDB_RESMAP_SECTION;
  while ( (ptr = infile.NextLine()) != 0 ) {
    ArgList argline( ptr, " " );
    // Check for section change first
    if (argline.Nargs() < 1) {
      if (section == PDB_RESMAP_SECTION) {
        //mprintf("DEBUG: Section change.\n");
        section = PDB_ATOMMAP_SECTION;
      } else if (section == PDB_ATOMMAP_SECTION) {
        section = PDB_LINKAGE_RES_SECTION;
      }
    } else if (argline[0][0] != '#') {
      // Skipping comments, read sections
      if (section == PDB_RESMAP_SECTION) {
        // OLD: "<Name>" <glycam reschar> <pdb resname list>
        // NEW: <res> <glycam> <form> <chirality> <ring type> <full name>
        SugarToken sToken;
        std::string sResName = sToken.SetFromLine(argline);
        if (sResName.empty()) {
          mprinterr("Error: Could not parse residue map section of '%s'\n",
                    infile.Filename().full());
          return 1;
        }
        std::pair<MapType::iterator,bool> ret =
          pdb_to_glycam_.insert( PairType(sResName, sToken) );
        if (!ret.second) {
          mprinterr("Error: Duplicate residue name '%s' in residue map section of '%s'\n",
                    sResName.c_str(), infile.Filename().full());
          return 1;
        }
      } else if (section == PDB_ATOMMAP_SECTION) {
        // <glycam reschar list> <PDB atomname to glycam atomname pair> ...
        if (argline.Nargs() < 2) {
          mprinterr("Error: Expected at least 2 columns in '%s' atom map section, got %i\n",
                    infile.Filename().full(), argline.Nargs());
          mprinterr("Error: %s\n", ptr);
          return 1;
        }
        // TODO handle glycam res names with > 1 char
        ArgList glycamnames( argline[0], "," );
        if (glycamnames.Nargs() < 1) {
          mprinterr("Error: No Glycam names found.\n");
          mprinterr("Error: %s\n", ptr);
          return 1;
        }
        int glycam_map_idx = (int)pdb_glycam_name_maps_.size();
        pdb_glycam_name_maps_.push_back(NameMapType());
        pdb_glycam_name_maps_A_.push_back(NameMapType());
        pdb_glycam_name_maps_B_.push_back(NameMapType());
        NameMapType& currentMap  = pdb_glycam_name_maps_.back();
        NameMapType& currentMapA = pdb_glycam_name_maps_A_.back();
        NameMapType& currentMapB = pdb_glycam_name_maps_B_.back();
        for (int col = 1; col < argline.Nargs(); col++) {
          ArgList namepair( argline[col], "," );
          NameMapType* currentMapPtr = &currentMap;
          if (namepair.Nargs() == 3) {
            // This name mapping is for a particular anomeric form
            if (namepair[2] == "A")
              currentMapPtr = &currentMapA;
            else if (namepair[2] == "B")
              currentMapPtr = &currentMapB;
            else {
              mprinterr("Error: For name pair, third arg should only be A or B: %s\n", ptr);
              return 1;
            }
          } else if (namepair.Nargs() != 2) {
            mprinterr("Error: Expected only 2 names for name pair, got %i\n", namepair.Nargs());
            mprinterr("Error: %s\n", ptr);
            return 1;
          }
          currentMapPtr->insert( NamePairType(NameType(namepair[0]), NameType(namepair[1])) );
        } // END loop over name pair columns
        // Map will be for each glycam res
        for (ArgList::const_iterator gres = glycamnames.begin(); gres != glycamnames.end(); ++gres)
          glycam_res_idx_map_.insert( ResIdxPairType( *gres, glycam_map_idx ) );
      } else if (section == PDB_LINKAGE_RES_SECTION) {
        // <pdb linkage res name> <glycam linkage res name>
        if (argline.Nargs() != 2) {
          mprinterr("Error: Expected only 2 columns in '%s' linkage res map section, got %i\n",
                    infile.Filename().full(), argline.Nargs());
          mprinterr("Error: %s\n", ptr);
        }
        pdb_glycam_linkageRes_map_.insert( NamePairType(NameType(argline[0]),
                                                        NameType(argline[1])) );
      }
    } // END not comment
  } // END loop over file
  infile.CloseFile();

  return 0;
}

// -----------------------------------------------------------------------------
