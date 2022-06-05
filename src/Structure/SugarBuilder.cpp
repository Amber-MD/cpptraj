#include "SugarBuilder.h"
#include "../ArgList.h"
#include "../CpptrajFile.h"
#include "../CpptrajStdio.h"
#include "../TorsionRoutines.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
SugarBuilder::SugarBuilder(int debugIn) : debug_(debugIn) {}

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
/// \return True if the given tgt atom is in the given array
static inline bool AtomIsInArray(std::vector<int> const& RingAtoms, int tgt)
{
  for (std::vector<int>::const_iterator it = RingAtoms.begin(); it != RingAtoms.end(); ++it)
    if (*it == tgt) return true;
  return false;
}

/// Recursive function for finding and recording all carbons
static void Find_Carbons(int atm, Topology const& topIn, std::vector<bool>& Visited,
                         std::vector<int>& remainingChainCarbons)
{
  remainingChainCarbons.push_back( atm );
  Visited[atm] = true;
  // Follow all carbons bonded to this atom
  for (Atom::bond_iterator bat = topIn[atm].bondbegin(); bat != topIn[atm].bondend(); ++bat)
  {
    if (topIn[*bat].Element() == Atom::CARBON && !Visited[*bat]) {
      Find_Carbons( *bat, topIn, Visited, remainingChainCarbons );
    }
  }
}

/** Find remaining non-ring carbons in chain starting from ring end atom. */
int SugarBuilder::FindRemainingChainCarbons(Iarray& remainingChainCarbons,
                                                   int start_c, Topology const& topIn, int rnum,
                                                   Iarray const& RingAtoms)
const
{
  Residue const& res = topIn.Res(rnum);
  std::vector<bool> Visited(topIn.Natom(), true);
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
    if (!AtomIsInArray(RingAtoms, at))
      Visited[at] = false;

  for (Atom::bond_iterator bat = topIn[start_c].bondbegin();
                           bat != topIn[start_c].bondend();
                         ++bat)
  {
    if ( !Visited[*bat] && topIn[*bat].Element() == Atom::CARBON )
      Find_Carbons(*bat, topIn, Visited, remainingChainCarbons);
  }
  return 0;
}
// -----------------------------------------------------------------------------
// Torsion routines
/// \return Position of given tgt atom in the array (if it is in the given array)
static inline int AtomIdxInArray(std::vector<int> const& ChainAtoms, int tgt)
{
  for (std::vector<int>::const_iterator it = ChainAtoms.begin(); it != ChainAtoms.end(); ++it)
    if (*it == tgt)
      return (int)(it - ChainAtoms.begin());
  return -1;
}

/** Determine torsion around the anomeric carbon. */
int SugarBuilder::CalcAnomericTorsion(double& torsion,
                                             int anomeric_atom, int ring_oxygen_atom,
                                             int rnum,
                                             Iarray const& RingAtoms,
                                             Topology const& topIn, Frame const& frameIn)
const
{
  if (debug_ > 0) {
    mprintf("\t  Anomeric carbon             : %s\n", topIn.ResNameNumAtomNameNum(anomeric_atom).c_str());
    mprintf("\t  Ring oxygen atom            : %s\n", topIn.ResNameNumAtomNameNum(ring_oxygen_atom).c_str());
  }
  int anomeric_atom_X = -1;
  int anomeric_atom_C = -1;
  // By definition the anomeric atom should be the first ring atom TODO catch size==1?
  anomeric_atom_C = RingAtoms[1];
  if (anomeric_atom_C == -1) {
    mprinterr("Error: Next ring atom after anomeric C could not be identified.\n");
    return 1;
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric C ring substituent : %s\n",
            topIn.ResNameNumAtomNameNum(anomeric_atom_C).c_str());
  // Get the substituent of the anomeric C (e.g. C1) that is a non-ring atom, non hydrogen 
  for ( Atom::bond_iterator bat = topIn[anomeric_atom].bondbegin();
                            bat != topIn[anomeric_atom].bondend();
                          ++bat )
  {
    if ( *bat != ring_oxygen_atom &&
         topIn[*bat].Element() != Atom::HYDROGEN &&
         !AtomIsInArray(RingAtoms, *bat) )
    {
      if (anomeric_atom_X != -1) {
        // If there are two non-ring, non-hydrogen substituents, prioritize
        // the one that is part of this residue.
        bool bat_in_res = (topIn[*bat].ResNum() == rnum);
        bool X_in_res   = (topIn[anomeric_atom_X].ResNum() == rnum);
        if ( (bat_in_res && X_in_res) || (!bat_in_res && !X_in_res) ) {
          // Both in or both out of residue. Choose based on atomic number.
          if (topIn[*bat].AtomicNumber() == topIn[anomeric_atom_X].AtomicNumber()) {
            mprinterr("Error: Two potential substituents for anomeric carbon: %s and %s\n",
                      topIn.ResNameNumAtomNameNum(*bat).c_str(),
                      topIn.ResNameNumAtomNameNum(anomeric_atom_X).c_str());
            return 1;
          } else if (topIn[*bat].AtomicNumber() > topIn[anomeric_atom_X].AtomicNumber()) {
            anomeric_atom_X = *bat;
          }
        } else if (bat_in_res) {
          anomeric_atom_X = *bat;
        }
      } else
        anomeric_atom_X = *bat;
    }
  }

  if (anomeric_atom_X == -1) {
    // If the Cx (C1 substituent, usually a different residue) index is
    // not found this usually means missing inter-residue bond.
    // Alternatively, this could be an isolated sugar missing an -OH
    // group, so make this non-fatal.
    mprintf("Warning: Anomeric C non-ring substituent could not be identified.\n");
    //        "Warning: This can happen if the sugar is bonded to something that\n"
    //        "Warning:  is missing, e.g. a -OH group. In that case the coordinates\n"
    //        "Warning   for the missing atoms may need to be generated.\n");
    return -1;
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric X substituent      : %s\n",
            topIn.ResNameNumAtomNameNum(anomeric_atom_X).c_str());

  torsion = Torsion( frameIn.XYZ(ring_oxygen_atom), frameIn.XYZ(anomeric_atom),
                     frameIn.XYZ(anomeric_atom_C), frameIn.XYZ(anomeric_atom_X) );
  if (debug_ > 0)
    mprintf("DEBUG: Anomeric torsion %s-%s-%s-%s= %f\n",
            *(topIn[ring_oxygen_atom].Name()),
            *(topIn[anomeric_atom].Name()),
            *(topIn[anomeric_atom_C].Name()),
            *(topIn[anomeric_atom_X].Name()),
            torsion * Constants::RADDEG);
  return 0;
}

/** Determine torsion around anomeric reference carbon. */
int SugarBuilder::CalcAnomericRefTorsion(double& torsion,
                                                int ano_ref_atom, int ring_oxygen_atom,
                                                int ring_end_atom, Iarray const& RingAtoms,
                                                Topology const& topIn, Frame const& frameIn)
const
{
  if (debug_ > 0)
    mprintf("\t  Anomeric ref carbon                   : %s\n",
            topIn.ResNameNumAtomNameNum(ano_ref_atom).c_str());
  //      ano_ref_atom_Y
  //           |
  //      ano_ref_atom
  //       |        |
  // ano_ref_atom_0 ano_ref_atom_1
  int ano_ref_atom_Y = -1;
  int ano_ref_atom_0 = -1;
  int ano_ref_atom_1 = -1;
  // This will be the index of the anomeric atom in the RingAtoms array
  int ar_index = -1;
  // Find ring atom that precedes the anomeric reference atom TODO catch size==1?
  for (unsigned int idx = 1; idx != RingAtoms.size(); idx++) {
    if (RingAtoms[idx] == ano_ref_atom) {
        ar_index = (int)idx;
        ano_ref_atom_0 = RingAtoms[idx-1];
        break;
    }
  }
  if (ano_ref_atom_0 == -1) {
    mprinterr("Error: Anomeric reference ring C previous ring atom could not be identified.\n");
    return 1;
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric reference previous ring atom : %s\n",
            topIn.ResNameNumAtomNameNum(ano_ref_atom_0).c_str());
  // If the anomeric reference atom is the ring end atom then ano_ref_atom_1
  // is the ring oxygen.
  if (ano_ref_atom == ring_end_atom) {
    ano_ref_atom_1 = ring_oxygen_atom;
  } else {
    // Anomeric reference atom is somewhere before the ring end atom.
    ano_ref_atom_1 = RingAtoms[ar_index+1];
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric reference next atom          : %s\n",
            topIn.ResNameNumAtomNameNum(ano_ref_atom_1).c_str());
  // Get non-hydrogen substituent of anomeric ref (e.g. C5) that 
  // is not part of the ring (e.g. C6).
  for ( Atom::bond_iterator bat = topIn[ano_ref_atom].bondbegin();
                            bat != topIn[ano_ref_atom].bondend();
                          ++bat )
  {
    if ( *bat != ring_oxygen_atom &&
         topIn[*bat].Element() != Atom::HYDROGEN &&
         !AtomIsInArray(RingAtoms, *bat) )
    {
      if (ano_ref_atom_Y != -1) {
        mprinterr("Error: Two potential non-ring substituents for anomeric ref: %s and %s\n",
                  topIn.ResNameNumAtomNameNum(*bat).c_str(),
                  topIn.ResNameNumAtomNameNum(ano_ref_atom_Y).c_str());
        return 1;
      }
      ano_ref_atom_Y = *bat;
    }
  }
  if (ano_ref_atom_Y == -1) {
    mprinterr("Error: Anomeric reference Y substituent could not be identified.\n");
    return 1;
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric reference substituent        : %s\n",
            topIn.ResNameNumAtomNameNum(ano_ref_atom_Y).c_str());

  torsion = Torsion( frameIn.XYZ(ano_ref_atom_0),   frameIn.XYZ(ano_ref_atom),
                     frameIn.XYZ(ano_ref_atom_1),   frameIn.XYZ(ano_ref_atom_Y) );
  if (debug_ > 0)
    mprintf("DEBUG: Anomeric reference torsion %s-%s-%s-%s= %f\n",
            *(topIn[ano_ref_atom_0].Name()),
            *(topIn[ano_ref_atom].Name()),
            *(topIn[ano_ref_atom_1].Name()),
            *(topIn[ano_ref_atom_Y].Name()),
            torsion * Constants::RADDEG);
  return 0;
}

/** Determine torsion around the configurational carbon. 
  * Calculate torsion around config. carbon C as:
  *   C0-C-Z-C1
  * where C0 is the carbon preceding C in the chain, C1 is the carbon
  * after C in the chain, and Z is the non-hydrogen substituent of C
  * with highest priority. Do it this way to be consistent with how
  * CalcAnomericRefTorsion orders the atoms.
  */
int SugarBuilder::CalcConfigCarbonTorsion(double& torsion, int config_carbon,
                                                 Iarray const& carbon_chain,
                                                 Topology const& topIn,
                                                 Frame const& frameIn)
const
{
  int atom_c0 = -1;
  int atom_c1 = -1;
  int atom_z  = -1;
  // Get c0 and c1
  int c_idx = AtomIdxInArray(carbon_chain, config_carbon);
  if (c_idx < 1) {
    mprinterr("Error: Could not determine carbon before config. C '%s'\n",
              topIn.ResNameNumAtomNameNum(config_carbon).c_str());
    return 1;
  }
  atom_c0 = carbon_chain[c_idx-1];
  if ((unsigned int)c_idx+1 >= carbon_chain.size()) {
    mprinterr("Error: Could not determine carbon after config. C '%s'\n",
              topIn.ResNameNumAtomNameNum(config_carbon).c_str());
    return 1;
  }
  atom_c1 = carbon_chain[c_idx+1];

  for (Atom::bond_iterator bat = topIn[config_carbon].bondbegin();
                           bat != topIn[config_carbon].bondend(); ++bat)
  {
    if (topIn[*bat].Element() != Atom::HYDROGEN &&
        !AtomIsInArray(carbon_chain, *bat))
    {
      if (atom_z == -1)
        atom_z = *bat;
      else if (topIn[*bat].AtomicNumber() > topIn[atom_z].AtomicNumber())
        atom_z = *bat;
    }
  }
  if (atom_z == -1) {
    mprinterr("Error: Could not determine substituent for config. C '%s'\n",
              topIn.ResNameNumAtomNameNum(config_carbon).c_str());
    return 1;
  }

  torsion = Torsion( frameIn.XYZ(atom_c0),
                     frameIn.XYZ(config_carbon),
                     frameIn.XYZ(atom_z),
                     frameIn.XYZ(atom_c1) );
  if (debug_ > 0)
    mprintf("DEBUG: Config. C torsion %s-%s-%s-%s= %f\n",
            *(topIn[atom_c0].Name()),
            *(topIn[config_carbon].Name()),
            *(topIn[atom_z].Name()),
            *(topIn[atom_c1].Name()),
            torsion*Constants::RADDEG);
  return 0;
}


