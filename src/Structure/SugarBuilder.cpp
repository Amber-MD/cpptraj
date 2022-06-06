#include "SugarBuilder.h"
#include "FxnGroupBuilder.h"
#include "ResStatArray.h"
#include "StructureRoutines.h"
#include "Sugar.h"
#include "SugarLinkAtoms.h"
#include "../ArgList.h"
#include "../CharMask.h" // for linkage determination
#include "../Chirality.h"
#include "../CpptrajFile.h"
#include "../CpptrajStdio.h"
#include "../DistRoutines.h"
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

// -----------------------------------------------------------------------------
/** \return true if residue is on the list of residues with missing heteroatoms. */
static bool res_missing_het(int rnum, Topology const& topIn) {
  Residue const& res = topIn.Res(rnum);
  //mprintf("DEBUG: TGT %s %i '%c' '%c'\n", *(res.Name()), res.OriginalResNum(), res.Icode(), res.ChainId());
  for (std::vector<Residue>::const_iterator het = topIn.MissingHet().begin();
                                            het != topIn.MissingHet().end(); ++het)
  {
    //mprintf("DEBUG: HET %s %i '%c' '%c'\n", *(het->Name()), het->OriginalResNum(), het->Icode(), het->ChainId());
    if ( het->OriginalResNum() == res.OriginalResNum() &&
         het->Icode()          == res.Icode() &&
         het->ChainId()        == res.ChainId() &&
         het->Name()           == res.Name() )
      return true;
  }
  return false;
}

/// Recursive function for following bonds of an atom to a target atom
static void FollowBonds(int atm, Topology const& topIn, int idx, std::vector<int>& ring_atoms, int tgt_atom, std::vector<bool>& Visited, bool& found)
{
  Visited[atm] = true;
  int rnum = topIn[atm].ResNum();
  //for (int i = 0; i != idx; i++) // DEBUG
  //  mprintf("\t"); // DBEUG
  //mprintf("At atom %s\n", topIn.ResNameNumAtomNameNum(atm).c_str()); // DEBUG
  ring_atoms[idx] = atm;
  // Assume we have started at the target atom
  if (idx > 0 && atm == tgt_atom) {
    found = true;
    return;
  }
  // Follow all atoms bonded to this atom
  for (Atom::bond_iterator bat = topIn[atm].bondbegin(); bat != topIn[atm].bondend(); ++bat)
  {
    if (topIn[*bat].ResNum() == rnum &&
        topIn[*bat].Element() == Atom::CARBON &&
        !Visited[*bat])
    {
      FollowBonds( *bat, topIn, idx+1, ring_atoms, tgt_atom, Visited, found );
      if (found) return;
    }
  }
}

/** Identify sugar oxygen, anomeric and ref carbons, and ring atoms. */
Sugar SugarBuilder::IdSugarRing(int rnum, Topology const& topIn)
const
{
  Residue const& res = topIn.Res(rnum);
  bool residue_missing_atoms = res_missing_het(rnum, topIn);
  if (residue_missing_atoms)
    mprintf("\tResidue '%s' is missing atoms.\n", topIn.TruncResNameOnumId(rnum).c_str());

  // Determine candidates for ring oxygen atoms. 
  Iarray potentialRingStartAtoms;
  // This array will hold the 2 flanking carbons for potential ring oxygen atoms.
  std::vector<std::pair<int,int>> ringAtomCarbons;
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
  {
    Atom const& currentAtom = topIn[at];
    // Try to identify the sugar ring oxygen. Candidate atoms are oxygens
    // bonded to two carbon atoms in the same residue.
    if (currentAtom.Element() == Atom::OXYGEN) {
      if (currentAtom.Nbonds() > 1) {
        Iarray c_atoms;
        for (Atom::bond_iterator bat = currentAtom.bondbegin();
                                 bat != currentAtom.bondend(); ++bat)
        {
          if (topIn[*bat].Element() == Atom::CARBON &&
              topIn[*bat].ResNum() == rnum)
            c_atoms.push_back( *bat );
        }
        if ( c_atoms.size() == 2 ) {
          potentialRingStartAtoms.push_back( at );
          ringAtomCarbons.push_back( std::pair<int,int>(c_atoms[0], c_atoms[1]) );
        }
      }
    }
  }

  if (potentialRingStartAtoms.empty()) {
    mprintf("Warning: Ring oxygen could not be identified for %s\n",
            topIn.TruncResNameOnumId(rnum).c_str());
    return Sugar(Sugar::MISSING_O, res.FirstAtom());
  }

  // Use the previously-set up AtomMap to help determine stereocenters
  std::vector<bool> atomIsChiral;
  atomIsChiral.reserve( res.NumAtoms() );
  // Since we cannot be certain there will be hydrogens, cannot rely
  // on the AtomMap chiral designations (which assumes hydrogens).
  // Make a chiral center a carbon with at least 3 bonds, all must be
  // to different kinds of atoms.
  int resat = res.FirstAtom();
  for (int iat = 0; iat != res.NumAtoms(); iat++, resat++)
  {
    bool chiral = false;
    if (topIn[resat].Element() == Atom::CARBON && topIn[resat].Nbonds() > 2) {
      if (debug_ > 0)
        mprintf("DEBUG: Atom '%s' potential chiral\n", topIn.TruncResAtomNameNum(resat).c_str());
      chiral = true;
      for (Atom::bond_iterator bat1 = topIn[resat].bondbegin();
                               bat1 != topIn[resat].bondend(); ++bat1)
      {
        std::string unique1 = myMap_[*bat1].Unique();
        for (Atom::bond_iterator bat2 = bat1 + 1; bat2 != topIn[resat].bondend(); ++bat2)
        {
          std::string unique2 = myMap_[*bat2].Unique();
          if (unique1 == unique2) {
            // At least two of the atoms bonded to this atom look the same. Not chiral.
            if (debug_ > 1)
              mprintf("DEBUG: unique strings match %s='%s' %s='%s'\n", *(topIn[*bat1].Name()), unique1.c_str(), *(topIn[*bat2].Name()), unique2.c_str());
            chiral = false;
            break;
          }
        } // END inner loop over bonded atoms
        if (!chiral) break;
      } // END outer loop over bonded atoms
    } // END atom is carbon with > 2 bonds
    
    atomIsChiral.push_back( chiral );
    if (debug_ > 0)
      mprintf("DEBUG: Atom '%s' isChiral= %i\n",
              topIn.TruncResAtomNameNum(resat).c_str(),
              (int)atomIsChiral.back());
  }

  // Using the potential ring start atoms, see if we can actually complete
  // a ring. Identify important atoms as well.
  // This will indicate ring direction
//  int ring_direction = 0;
  // Ring end atom is the last atom in the ring
  int ring_end_atom = -1;
  // The anomeric carbon is the carbon that was part of the carbonyl group
  // in the straight chain. It is therefore typically the carbon with fewer
  // bonds to other carbons.
  int anomeric_atom = -1;    // e.g. C1
  // The anomeric reference carbon is the stereocenter farthest from the
  // anomeric carbon in the ring.
  int ano_ref_atom = -1;     // e.g. C5
  // Ring oxygen atom
  int ring_oxygen_atom = -1; // e.g. O5
  // This will hold the index of the highest stereocenter, e.g. C5
  int highest_stereocenter = -1;
  // This will hold ring atoms, not including the ring oxygen.
  std::vector<Iarray> Ring_Atoms;
  // This will hold carbon chain atoms starting from the anomeric carbon
  Iarray carbon_chain;

  // Out of the potential ring start atoms, see which ones are actually
  // part of a ring. Potential ring start atoms only have 2 bonds,
  // each one to a carbon.
  int ring_atom_idx = -1;
  std::vector<std::pair<int,int>>::iterator catoms = ringAtomCarbons.begin();
  for (Iarray::const_iterator ringat = potentialRingStartAtoms.begin();
                              ringat != potentialRingStartAtoms.end();
                            ++ringat, ++catoms)
  {
    if (debug_ > 0)
      mprintf("DEBUG: Ring start '%s'\n", topIn.ResNameNumAtomNameNum(*ringat).c_str());
    // Mark all atoms as visited except this residue (minus the ring start).
    std::vector<bool> Visited( topIn.Natom(), true );
    for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
      if (at != *ringat)
        Visited[at] = false;
    Iarray ring_atoms( topIn.Res(rnum).NumAtoms(), -1 );
    // Will be set true if complete ring can be found
    bool ring_complete = false;

    // Since we have already established that *ringat is an oxygen bonded
    // to two carbons, just start at the first carbon to see if we can
    // get to the second carbon.
    int c_beg, c_end;
    if (catoms->first < catoms->second) {
      c_beg = catoms->first;
      c_end = catoms->second;
    } else {
      c_beg = catoms->second;
      c_end = catoms->first;
    }
    // Try to ascertain which carbon might be the anomeric carbon (i.e. the
    // carbon that originally started the chain). Tie goes to lower index.
    int c_beg_bonds_to_C = 0;
    for (Atom::bond_iterator bat = topIn[c_beg].bondbegin(); bat != topIn[c_beg].bondend(); ++bat)
      if (topIn[*bat].Element() == Atom::CARBON)
        c_beg_bonds_to_C++;
    int c_end_bonds_to_C = 0;
    for (Atom::bond_iterator bat = topIn[c_end].bondbegin(); bat != topIn[c_end].bondend(); ++bat)
      if (topIn[*bat].Element() == Atom::CARBON)
        c_end_bonds_to_C++;
    if (debug_ > 0)
      mprintf("DEBUG:\t(%s bonds to C= %i, %s bonds to C = %i)\n", // DEBUG
              topIn.ResNameNumAtomNameNum(c_beg).c_str(), c_beg_bonds_to_C,
              topIn.ResNameNumAtomNameNum(c_end).c_str(), c_end_bonds_to_C);
    if (c_beg_bonds_to_C <= c_end_bonds_to_C) {
      anomeric_atom = c_beg;
      ring_end_atom = c_end;
//      ring_direction = 1;
    } else {
      anomeric_atom = c_end;
      ring_end_atom = c_beg;
//      ring_direction = -1;
    }
//    mprintf("DEBUG: Potential Ring direction= %i\n", ring_direction);
    catoms->first = anomeric_atom;
    catoms->second = ring_end_atom;

    FollowBonds(anomeric_atom, topIn, 0, ring_atoms,
                ring_end_atom, Visited, ring_complete);
    if (debug_ > 0)
      mprintf("DEBUG: Potential ring start atom %s, Ring complete = %i",
              topIn.ResNameNumAtomNameNum(*ringat).c_str(), (int)ring_complete);
    // Create empty array for ring atoms
    Ring_Atoms.push_back(Iarray());
    if (ring_complete) {
      // Able to complete the cycle.
      if (ring_atom_idx == -1)
        ring_atom_idx = (int)(ringat - potentialRingStartAtoms.begin());
      else {
        mprinterr("Error: Multiple potential ring atoms: %s and %s\n",
                  topIn.ResNameNumAtomNameNum(potentialRingStartAtoms[ring_atom_idx]).c_str(),
                  topIn.ResNameNumAtomNameNum(*ringat).c_str());
        return Sugar(Sugar::MULTIPLE_O, res.FirstAtom());
      }
      // Place the ring atoms into an array without the terminating -1
      Iarray& RA = Ring_Atoms.back();
      if (debug_ > 0) mprintf(" :"); // DEBUG
      for (Iarray::const_iterator it = ring_atoms.begin(); it != ring_atoms.end(); ++it)
      {
        if (debug_ > 0) mprintf(" %i", *it + 1);
        if (*it == -1) break;
        RA.push_back( *it );
      }
    }
    if (debug_ > 0) mprintf("\n"); // DEBUG
  } // END loop over potential ring atoms
  if (ring_atom_idx == -1) {
    mprinterr("Error: Sugar ring oxygen could not be identified.\n");
    return Sugar(Sugar::MISSING_O, res.FirstAtom());
  }

  ring_oxygen_atom = potentialRingStartAtoms[ring_atom_idx];
  anomeric_atom = ringAtomCarbons[ring_atom_idx].first;
  ring_end_atom = ringAtomCarbons[ring_atom_idx].second;
  Iarray const& RA = Ring_Atoms[ring_atom_idx];

  // Find anomeric reference atom. Start at ring end and work down to anomeric atom
  for (Iarray::const_iterator arat = RA.end() - 1; arat != RA.begin(); --arat)
    if (atomIsChiral[*arat - topIn.Res(rnum).FirstAtom()]) {
      ano_ref_atom = *arat;
      break;
    }

  // Get complete chain starting from the anomeric carbon
  carbon_chain = RA;
  if (FindRemainingChainCarbons(carbon_chain, ring_end_atom, topIn, rnum, RA)) {
    mprinterr("Error: Could not find remaining chain carbons.\n");
    return Sugar(Sugar::MISSING_CHAIN, res.FirstAtom());
  }
  if (debug_ > 0) {
    mprintf("DEBUG: Complete carbon chain (from anomeric carbon):\n");
    for (Iarray::const_iterator it = carbon_chain.begin(); it != carbon_chain.end(); ++it)
      mprintf("\t\t%s\n", topIn.ResNameNumAtomNameNum(*it).c_str());
  }
  // See if there is chain prior to anomeric carbon
  Iarray previous_chain;
  if (FindRemainingChainCarbons(previous_chain, anomeric_atom, topIn, rnum, RA)) {
    mprinterr("Error: Could not find previous chain carbons.\n");
    return Sugar(Sugar::MISSING_CHAIN, res.FirstAtom());
  }
  if (!previous_chain.empty()) {
    if (debug_ > 0) {
      mprintf("DEBUG: Previous carbon chain (from anomeric carbon):\n");
      for (Iarray::const_iterator it = previous_chain.begin(); it != previous_chain.end(); ++it)
        mprintf("\t\t%s\n", topIn.ResNameNumAtomNameNum(*it).c_str());
    }
    for (Iarray::const_iterator it = carbon_chain.begin(); it != carbon_chain.end(); ++it)
      previous_chain.push_back( *it );
    carbon_chain = previous_chain;
    if (debug_ > 0) {
      mprintf("DEBUG: Complete carbon chain:\n");
      for (Iarray::const_iterator it = carbon_chain.begin(); it != carbon_chain.end(); ++it)
        mprintf("\t\t%s\n", topIn.ResNameNumAtomNameNum(*it).c_str());
    }
  } 
  // Get the index of the highest stereocenter
  for (Iarray::const_iterator it = carbon_chain.begin(); it != carbon_chain.end(); ++it)
  {
    if (atomIsChiral[*it - topIn.Res(rnum).FirstAtom()])
      highest_stereocenter = *it;
  }
  if (debug_ > 0)
    mprintf("DEBUG: Index of highest stereocenter: %s\n",
            topIn.ResNameNumAtomNameNum(highest_stereocenter).c_str());

  if (ano_ref_atom == -1) {
    mprinterr("Error: Anomeric reference atom could not be identified.\n");
    return Sugar(Sugar::MISSING_ANO_REF, ring_oxygen_atom, anomeric_atom, RA, carbon_chain);
  }
  if (highest_stereocenter == -1) {
    mprinterr("Error: Highest stereocenter atom could not be identified.\n");
    return Sugar(Sugar::MISSING_CONFIG, ring_oxygen_atom, anomeric_atom, RA, carbon_chain);
  }

//  if (!ring_complete || RA.empty() || ring_oxygen_atom == -1) {
//    mprinterr("Error: Sugar ring atoms could not be identified.\n");
//    stat = ID_ERR;
//    return Sugar(res.FirstAtom());
//  }
  if (debug_ > 0)
    mprintf("\t  Ring oxygen         : %s\n", topIn.ResNameNumAtomNameNum(ring_oxygen_atom).c_str());
  return Sugar(ring_oxygen_atom, anomeric_atom, ano_ref_atom, highest_stereocenter,
               RA, carbon_chain, residue_missing_atoms);
}

// -----------------------------------------------------------------------------
/** Change PDB atom names in residue to glycam ones. */
int SugarBuilder::ChangePdbAtomNamesToGlycam(std::string const& resCode, Residue const& res,
                                             Topology& topIn, SugarToken::FormTypeEnum form)
const
{
  // Get the appropriate map
  ResIdxMapType::const_iterator resIdxPair = glycam_res_idx_map_.find( resCode );
  if (resIdxPair == glycam_res_idx_map_.end()) {
    // No map needed for this residue
    //mprintf("DEBUG: No atom map for residue '%s'.\n", resCode.c_str());
    return 0;
  }
  NameMapType const& currentMap = pdb_glycam_name_maps_[resIdxPair->second];
  NameMapType const* currentMapAB;
  if (form == SugarToken::ALPHA)
    currentMapAB = &(pdb_glycam_name_maps_A_[resIdxPair->second]);
  else
    currentMapAB = &(pdb_glycam_name_maps_B_[resIdxPair->second]);
  // Change PDB names to Glycam ones
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
  {
    NameMapType::const_iterator namePair = currentMapAB->find( topIn[at].Name() );
    if (namePair != currentMapAB->end())
      ChangeAtomName( topIn.SetAtom(at), namePair->second );
    else {
      namePair = currentMap.find( topIn[at].Name() );
      if (namePair != currentMap.end())
        ChangeAtomName( topIn.SetAtom(at), namePair->second );
    }
  }
  return 0;
}

// -----------------------------------------------------------------------------
/** Determine if anomeric carbon of furanose is up or down. */
int SugarBuilder::DetermineUpOrDown(SugarToken& stoken,
                                         Sugar const& sugar,
                                         Topology const& topIn, Frame const& frameIn)
const
{
  int cdebug;
  if (debug_ > 1)
    cdebug = 1;
  else
    cdebug = 0;
  Cpptraj::Chirality::ChiralType ctypeR = Cpptraj::Chirality::
                                          DetermineChirality(sugar.HighestStereocenter(),
                                                             topIn, frameIn, cdebug);
  if (ctypeR == Cpptraj::Chirality::ERR) {
    mprinterr("Error: Could not determine configuration for furanose.\n"); // TODO warn?
    return 1;
  }
  if (ctypeR == Cpptraj::Chirality::IS_R)
    stoken.SetChirality(SugarToken::IS_D);
  else
    stoken.SetChirality(SugarToken::IS_L);

  Cpptraj::Chirality::ChiralType ctypeA = Cpptraj::Chirality::
                                          DetermineChirality(sugar.AnomericAtom(),
                                                             topIn, frameIn, cdebug);
  if (ctypeA == Cpptraj::Chirality::ERR) {
    mprinterr("Error: Could not determine chirality around anomeric atom for furanose.\n"); // TODO warn?
    return 1;
  }

  if (ctypeR == ctypeA) {
    // Up, beta
    stoken.SetForm(SugarToken::BETA);
  } else {
    // Down, alpha
    stoken.SetForm(SugarToken::ALPHA);
  }
  return 0;
}

/** Determine anomeric form of the sugar. */
int SugarBuilder::DetermineAnomericForm(SugarToken& stoken,
                                             Sugar& sugarIn,
                                             Topology const& topIn, Frame const& frameIn)
const
{
  Sugar const& sugar = sugarIn;
  // For determining orientation around anomeric carbon need ring
  // oxygen atom and next carbon in the ring.
  double t_an;
//  ChiralRetType ac_chirality = CalcChiralAtomTorsion(t_an, anomeric_atom, topIn, frameIn);
//  mprintf("DEBUG: Based on t_an %s chirality is %s\n",
//          topIn.TruncResNameOnumId(rnum).c_str(), chiralStr[ac_chirality]);
  int ret = CalcAnomericTorsion(t_an, sugar.AnomericAtom(), sugar.RingOxygenAtom(),
                                sugar.ResNum(topIn),
                                sugar.RingAtoms(), topIn, frameIn);
  if (ret < 0) {
    // This means C1 X substituent missing; non-fatal.
    sugarIn.SetStatus( Sugar::MISSING_C1X );
    return 1; // TODO return 0?
  } else if (ret > 0) {
    // Error
    return 1; 
  }
  bool t_an_up = (t_an > 0);

  // For determining orientation around anomeric reference carbon need
  // previous carbon in the chain and either next carbon or ring oxygen.
  double t_ar;
//    ChiralRetType ar_chirality = CalcChiralAtomTorsion(t_ar, ano_ref_atom, topIn, frameIn);
//    mprintf("DEBUG: Based on t_ar %s chirality is %s\n",
//            topIn.TruncResNameOnumId(rnum).c_str(), chiralStr[ar_chirality]);
  if (CalcAnomericRefTorsion(t_ar, sugar.AnomericRefAtom(), sugar.RingOxygenAtom(), sugar.RingEndAtom(),
                             sugar.RingAtoms(), topIn, frameIn))
  {
    return 1; 
  }
  bool t_ar_up = (t_ar > 0);

  // If config. C is not the anomeric reference, need the previous
  // carbon in the chain, next carbon in the chain, and config. C
  // substituent.
  double t_cc;
  if (sugar.AnomericRefAtom() != sugar.HighestStereocenter()) {
    if (CalcConfigCarbonTorsion(t_cc, sugar.HighestStereocenter(),
                                sugar.ChainAtoms(), topIn, frameIn))
      return 1;
  } else
    t_cc = t_ar;
  bool t_cc_up = (t_cc > 0);

  // Determine index of anomeric atom (typically index 0 but not always).
  int aa_idx =  AtomIdxInArray(sugar.ChainAtoms(), sugar.AnomericAtom());
  int aa_pos = (aa_idx % 2);
  // Determine index of the anomeric reference atom in the chain.
  int ar_idx = AtomIdxInArray(sugar.ChainAtoms(), sugar.AnomericRefAtom());
  int cc_idx = AtomIdxInArray(sugar.ChainAtoms(), sugar.HighestStereocenter());

  // Determine form and chirality.
  // May need to adjust definitions based on the positions of the anomeric
  // reference and config. atoms in the sequence, which alternates.
  if ((ar_idx % 2) != aa_pos)
    t_ar_up = !t_ar_up;
  if ((cc_idx % 2) != aa_pos)
    t_cc_up = !t_cc_up;

  if ( debug_ > 0) {
    mprintf("DEBUG: Index of the anomeric reference atom is %i\n", ar_idx);
    mprintf("DEBUG: Index of the config. carbon atom is %i\n", cc_idx);
    mprintf("DEBUG: t_an_up=%i  t_ar_up=%i  t_cc_up=%i\n",
            (int)t_an_up, (int)t_ar_up, (int)t_cc_up);
  }

  // Same side is beta, opposite is alpha.
  if (t_an_up == t_ar_up) {
    stoken.SetForm(SugarToken::BETA); //form = IS_BETA;
    //mprintf("DEBUG: Form is Beta\n");
  } else {
    stoken.SetForm(SugarToken::ALPHA); //form = IS_ALPHA;
    //mprintf("DEBUG: Form is Alpha\n");
  }

  // By the atom ordering used by CalcAnomericRefTorsion and
  // CalcConfigCarbonTorsion, D is a negative (down) torsion.
  if (!t_cc_up)
    stoken.SetChirality(SugarToken::IS_D);
  else
    stoken.SetChirality(SugarToken::IS_L);

  return 0;
}

// -----------------------------------------------------------------------------

/** Determine linkages for the sugar. Non-sugar residues linked to sugars
  * with a recognized linkage will be marked as valid.
  */
std::string SugarBuilder::DetermineSugarLinkages(Sugar const& sugar, CharMask const& cmask,
                                                 Topology& topIn, ResStatArray& resStatIn,
//                                                        CpptrajFile* outfile,
                                                        std::set<BondType>& sugarBondsToRemove)
const
{
  int rnum = sugar.ResNum(topIn);
  Residue const& res = topIn.Res(rnum);

  // If we haven't identified carbon chain atoms we can't determine their
  // position and therefore the link code.
  if (sugar.ChainAtoms().empty()) {
    mprinterr("Error: No chain atoms determined for '%s', cannot determine link positions.\n",
              topIn.TruncResNameOnumId(rnum).c_str());
    resStatIn[rnum] = ResStatArray::SUGAR_NO_CHAIN_FOR_LINK;
    return std::string("");
  }

  // Use set to store link atoms so it is ordered by position
  SugarLinkAtoms linkages(debug_);

  // Bonds to non sugars to be removed since these will confuse tleap
  //BondArray bondsToRemove;

  // Loop over sugar atoms
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
  {
    // Ignore the sugar oxygen and any hydrogen atoms
    if (at == sugar.RingOxygenAtom() ||
        topIn[at].Element() == Atom::HYDROGEN)
      continue;
    // This will be the index of the carbon that is atom is bonded to (or of this atom itself).
    int atomChainPosition = -1;
    // Find position in carbon chain
    if (topIn[at].Element() == Atom::CARBON)
      atomChainPosition = AtomIdxInArray(sugar.ChainAtoms(), at);
    if (atomChainPosition == -1) {
      // Check if any bonded atoms are in carbon chain
      for (Atom::bond_iterator bat = topIn[at].bondbegin();
                               bat != topIn[at].bondend(); ++bat)
      {
        if (topIn[*bat].Element() == Atom::CARBON) {
          atomChainPosition = AtomIdxInArray(sugar.ChainAtoms(), *bat);
          if (atomChainPosition != -1) break;
        }
      }
    }

    // Check for bonds to other residues
    for (Atom::bond_iterator bat = topIn[at].bondbegin();
                             bat != topIn[at].bondend(); ++bat)
    {
      if (topIn[*bat].ResNum() != rnum) {
        // This atom is bonded to another residue.
        linkages.AddLinkAtom(at, atomChainPosition+1);
        // Check if the other residue is a sugar or not
        if (!cmask.AtomInCharMask(*bat)) {
          // Atom is bonded to non-sugar residue.
          mprintf("\t  Sugar %s %s (%i) bonded to non-sugar %s %s (%i) at position %i\n",
                  topIn.TruncResNameOnumId(rnum).c_str(), *(topIn[at].Name()), rnum+1,
                  topIn.TruncResNameOnumId(topIn[*bat].ResNum()).c_str(),
                  *(topIn[*bat].Name()), topIn[*bat].ResNum()+1, atomChainPosition+1);
          sugarBondsToRemove.insert( BondType(at, *bat, -1) );
          // Check if this is a recognized linkage to non-sugar
          Residue& pres = topIn.SetRes( topIn[*bat].ResNum() );
          NameMapType::const_iterator lname = pdb_glycam_linkageRes_map_.find( pres.Name() );
          if (lname != pdb_glycam_linkageRes_map_.end()) {
            if (debug_ > 0)
              mprintf("DEBUG: Link residue name for %s found: %s\n", *(lname->first), *(lname->second));
            ChangeResName( pres, lname->second );
            resStatIn[topIn[*bat].ResNum()] = ResStatArray::VALIDATED;
          
          } else if (pres.Name() == "ROH") {
            if (debug_ > 0)
              mprintf("DEBUG: '%s' is terminal hydroxyl.\n", *(pres.Name()));
            resStatIn[topIn[*bat].ResNum()] = ResStatArray::VALIDATED;
          } else if (pres.Name() == "SO3") {
            if (debug_ > 0)
              mprintf("DEBUG: '%s' is a sulfate group.\n", *(pres.Name()));
            resStatIn[topIn[*bat].ResNum()] = ResStatArray::VALIDATED;
          } else if (pres.Name() == "MEX") {
            if (debug_ > 0)
              mprintf("DEBUG: '%s' is a methyl group.\n", *(pres.Name()));
            resStatIn[topIn[*bat].ResNum()] = ResStatArray::VALIDATED;
          } else if (pres.Name() == "OME") {
            if (debug_ > 0)
              mprintf("DEBUG: '%s' is an O-methyl group.\n", *(pres.Name()));
            resStatIn[topIn[*bat].ResNum()] = ResStatArray::VALIDATED;
          } else if (hasGlycam_) {
            // Check if pres.Name() is already changed to linkage name
            for (NameMapType::const_iterator rname = pdb_glycam_linkageRes_map_.begin();
                                             rname != pdb_glycam_linkageRes_map_.end(); ++rname)
            {
              if (pres.Name() == rname->second) {
                if (debug_ > 0)
                  mprintf("DEBUG: Link residue for %s (%s) is already %s\n",
                          topIn.TruncResNameOnumId(topIn[*bat].ResNum()).c_str(),
                          *(rname->first), *(rname->second));
                resStatIn[topIn[*bat].ResNum()] = ResStatArray::VALIDATED;
                break;
              }
            }
          }
          if (resStatIn[topIn[*bat].ResNum()] != ResStatArray::VALIDATED) {
            mprintf("Warning: Unrecognized link residue %s, not modifying name.\n", *pres.Name());
            resStatIn[topIn[*bat].ResNum()] = ResStatArray::SUGAR_UNRECOGNIZED_LINK_RES;
          }
        } else {
          // Atom is bonded to sugar residue
          mprintf("\t  Sugar %s %s bonded to sugar %s %s at position %i\n",
                  topIn.TruncResNameOnumId(rnum).c_str(), *(topIn[at].Name()),
                  topIn.TruncResNameOnumId(topIn[*bat].ResNum()).c_str(),
                  *(topIn[*bat].Name()), atomChainPosition+1);
          // Also remove inter-sugar bonds since leap cant handle branching
          if (at < *bat)
            sugarBondsToRemove.insert( BondType(at, *bat, -1) );
          else
            sugarBondsToRemove.insert( BondType(*bat, at, -1) );
        }
      }
    } // END loop over bonded atoms

  } // END loop over sugar atoms

  // Determine linkage
  //if (debug_ > 0) {
    mprintf("\t  Link atoms:");
    for (SugarLinkAtoms::const_iterator it = linkages.begin();
                                        it != linkages.end(); ++it)
      mprintf(" %s(%i)", *(topIn[it->Idx()].Name()), it->Position());
    mprintf("\n");
  //}
  std::string linkcode;
  if (linkages.NoLinkAtoms()) {
    mprintf("\t  No linkages (may be missing atoms).\n");
    resStatIn[rnum] = ResStatArray::SUGAR_NO_LINKAGE;
    return linkcode;
  } else {
    linkcode = linkages.GlycamLinkageCode(topIn);
    if (debug_ > 0) mprintf("\t  Linkage code: %s\n", linkcode.c_str());
    if (linkcode.empty()) {
      mprinterr("Error: Unrecognized sugar linkage.\n");
      resStatIn[rnum] = ResStatArray::SUGAR_UNRECOGNIZED_LINKAGE;
      return linkcode;
    }
  }
  // Remove bonds to other residues 
  /*for (BondArray::const_iterator bnd = bondsToRemove.begin();
                                 bnd != bondsToRemove.end(); ++bnd)
  {
    LeapBond(bnd->A1(), bnd->A2(), topIn, outfile);
    topIn.RemoveBond(bnd->A1(), bnd->A2());
  }*/

  return linkcode;
}

/** Create a residue mask string for selecting Glycam-named sugar residues. */
std::string SugarBuilder::GenGlycamResMaskString() const {
  // Linkage codes
  // TODO this is a bit of a hack since not all linkage codes are valid for each sugar.
  // Should tighten this up later...
  typedef std::vector<std::string> Sarray;
  std::vector< std::string > linkageCodes;
  linkageCodes.push_back("0");
  linkageCodes.push_back("1");
  linkageCodes.push_back("2");
  linkageCodes.push_back("3");
  linkageCodes.push_back("4");
  linkageCodes.push_back("5");
  linkageCodes.push_back("6");
  linkageCodes.push_back("Z");
  linkageCodes.push_back("Y");
  linkageCodes.push_back("X");
  linkageCodes.push_back("W");
  linkageCodes.push_back("V");
  linkageCodes.push_back("U");
  linkageCodes.push_back("T");
  linkageCodes.push_back("S");
  linkageCodes.push_back("R");
  linkageCodes.push_back("Q");
  linkageCodes.push_back("P");
  // Loop over names
  std::string maskString;
  std::set< std::string > glycamResNames;
  for (MapType::const_iterator mit = pdb_to_glycam_.begin(); mit != pdb_to_glycam_.end(); ++mit)
  {
    SugarToken const& tkn = mit->second;
    std::pair<std::set<std::string>::iterator, bool> ret = glycamResNames.insert( tkn.GlycamCode() );
    if (ret.second == false) {
      if (debug_ > 1)
        mprintf("DEBUG: Already seen '%s', skipping.\n", tkn.GlycamCode().c_str());
      continue;
    }
    // L/D forms
    std::string glycamCodes[2];
    // Alpha/beta forms
    std::string formCodes[2];
    // D form is uppercase, which is the default.
    glycamCodes[0] = tkn.GlycamCode();
    // L form has lowercase
    for (std::string::const_iterator it = glycamCodes[0].begin(); it != glycamCodes[0].end(); ++it)
      glycamCodes[1] += tolower( *it );
    // Alpha/beta based on ring type
    if (tkn.RingType() == SugarToken::PYRANOSE) {
      formCodes[0] = "A";
      formCodes[1] = "B";
    } else if (tkn.RingType() == SugarToken::FURANOSE) {
      formCodes[0] = "D";
      formCodes[1] = "U";
    } else {
      mprinterr("Internal Error: Unhandled ring type in GenGlycamResMaskString().\n");
      return std::string("");
    }
    // Linkage codes
    // Create strings
    for (int ii = 0; ii != 2; ii++) {
      for (int jj = 0; jj != 2; jj++) {
        for (Sarray::const_iterator it = linkageCodes.begin(); it != linkageCodes.end(); ++it) {
          if (maskString.empty())
            maskString.assign( ":" + *it + glycamCodes[ii] + formCodes[jj] );
          else
            maskString.append( "," +  *it + glycamCodes[ii] + formCodes[jj] );
        }
      }
    }
  }

  return maskString;
}
// -----------------------------------------------------------------------------
/** Attempt to identify sugar residue, form, and linkages. */
int SugarBuilder::IdentifySugar(Sugar& sugarIn, Topology& topIn,
                                       Frame const& frameIn, CharMask const& cmask,
                                       ResStatArray& resStatIn,
                                       std::set<BondType>& sugarBondsToRemove)
{
  Sugar const& sugar = sugarIn;
  const std::string sugarName = topIn.TruncResNameOnumId(sugar.ResNum(topIn));
  const std::string resName = topIn.Res(sugar.ResNum(topIn)).Name().Truncated();
  //if (sugar.NotSet()) {
  //  mprintf("Warning: Sugar %s is not set up. Skipping sugar identification.\n", sugarName.c_str());
  //  return 0; // TODO return 1?
  //}

  int rnum = sugar.ResNum(topIn);
  Residue& res = topIn.SetRes(rnum);
  // Try to ID the base sugar type from the input name.
  //char resChar = ' ';

  MapType::const_iterator pdb_glycam;
  if (!hasGlycam_) {
    pdb_glycam = pdb_to_glycam_.find( res.Name() );
    //resChar = pdb_glycam->second;
  } else {
    mprintf("\tGlycam res: '%s'\n", resName.c_str());
    // Assume residue name is already glycam. Find it. Get sugar from 2nd char.
    std::string rChar(1, resName[1]);
    // Determine chirality from upper/lowercase
    SugarToken::ChirTypeEnum cType = SugarToken::UNKNOWN_CHIR;
    if (islower(rChar.front()))
      cType = SugarToken::IS_L;
    else
      cType = SugarToken::IS_D;
    // Now force uppercase
    rChar.assign( 1, toupper(rChar.front()) );
    std::string fChar(1, resName[2]);
    if (debug_ > 0)
      mprintf("DEBUG: Searching for glycam char '%s' '%s'\n", rChar.c_str(), fChar.c_str());
    // Get form from 3rd char
    SugarToken::FormTypeEnum fType = SugarToken::UNKNOWN_FORM;
    if (fChar == "U" || fChar == "B")
      fType = SugarToken::BETA;
    else if (fChar == "D" || fChar == "A")
      fType = SugarToken::ALPHA;
    for (pdb_glycam = pdb_to_glycam_.begin(); pdb_glycam != pdb_to_glycam_.end(); ++pdb_glycam)
      if (pdb_glycam->second.GlycamCode() == rChar &&
          pdb_glycam->second.Form() == fType &&
          pdb_glycam->second.Chirality() == cType)
        break;
  }

  if ( pdb_glycam == pdb_to_glycam_.end() ) { 
    mprinterr("Error: Could not identify sugar from residue name '%s'\n", *res.Name());
    return 1;
  }

  mprintf("\tSugar %s %s\n", sugarName.c_str(), pdb_glycam->second.InfoStr().c_str());

  SugarToken sugarInfo;
  int detect_err = 1;

  if (!sugar.NotSet()) {
    sugarInfo = SugarToken(sugar.RingType());

    // Determine alpha or beta and D or L
    if (sugarInfo.RingType() == SugarToken::FURANOSE)
      detect_err = DetermineUpOrDown(sugarInfo, sugar, topIn, frameIn);
    else if (sugarInfo.RingType() == SugarToken::PYRANOSE)
      detect_err = DetermineAnomericForm(sugarInfo, sugarIn, topIn, frameIn);
    else
      detect_err = 1;
  }

  if (detect_err != 0) {
    // This means there was an issue determining ring type, form,
    // chirality, or both.
    // Try to fill in info based on the residue name (pdb_glycam->second).
    if (sugarInfo.Form() == SugarToken::UNKNOWN_FORM) {
      mprintf("Warning: Could not determine anomer type from coordinates.\n");
      if (pdb_glycam->second.Form() == SugarToken::UNKNOWN_FORM)
        return 0;
      mprintf("Warning: Setting anomer type based on residue name.\n");
      sugarInfo.SetForm( pdb_glycam->second.Form() );
    }
    if (sugarInfo.Chirality() == SugarToken::UNKNOWN_CHIR) {
      mprintf("Warning: Could not determine configuration from coordinates.\n");
      if (pdb_glycam->second.Chirality() == SugarToken::UNKNOWN_CHIR)
        return 0;
      mprintf("Warning: Setting configuration based on residue name.\n");
      sugarInfo.SetChirality( pdb_glycam->second.Chirality() );
    }
    if (sugarInfo.RingType() == SugarToken::UNKNOWN_RING) {
      mprintf("Warning: Could not determine ring type from coordinates.\n");
      if (pdb_glycam->second.RingType() == SugarToken::UNKNOWN_RING)
        return 0;
      mprintf("Warning: Setting ring type based on residue name.\n");
      sugarInfo.SetRingType( pdb_glycam->second.RingType() );
    }
  }
  // Warn about form/chirality mismatches. If a mismatch occurs and the sugar
  // did not have all atoms present, defer to the name.
  if (pdb_glycam->second.Form() != SugarToken::UNKNOWN_FORM &&
      pdb_glycam->second.Form() != sugarInfo.Form())
  {
    mprintf("Warning: '%s' detected anomer type is %s but anomer type based on name is %s.\n",
            sugarName.c_str(), sugarInfo.formStr(),
            SugarToken::formStr(pdb_glycam->second.Form()) );
    if (sugarInfo.Form() != SugarToken::UNKNOWN_FORM)
      resStatIn[rnum] = ResStatArray::SUGAR_NAME_MISMATCH;
    if (useSugarName_) {
      mprintf("\tSetting anomer type based on residue name.\n");
      sugarInfo.SetForm( pdb_glycam->second.Form() );
    } else if (sugar.IsMissingAtoms()) {
      mprintf("Warning: Residue was missing atoms; setting anomer type based on residue name.\n");
      sugarInfo.SetForm( pdb_glycam->second.Form() );
    }
  }
  if (pdb_glycam->second.Chirality() != SugarToken::UNKNOWN_CHIR &&
      pdb_glycam->second.Chirality() != sugarInfo.Chirality())
  {
    mprintf("Warning: '%s' detected configuration is %s but configuration based on name is %s.\n",
            sugarName.c_str(), sugarInfo.chirStr(),
            SugarToken::chirStr(pdb_glycam->second.Chirality()) );
    if (sugarInfo.Chirality() != SugarToken::UNKNOWN_CHIR)
      resStatIn[rnum] = ResStatArray::SUGAR_NAME_MISMATCH;
    if (useSugarName_) {
      mprintf("\tSetting configuration based on residue name.\n");
      sugarInfo.SetChirality( pdb_glycam->second.Chirality() );
    } else if (sugar.IsMissingAtoms()) {
      mprintf("Warning: Residue was missing atoms; setting configuration based on residue name.\n");
      sugarInfo.SetChirality( pdb_glycam->second.Chirality() );
    }
  }
  if (pdb_glycam->second.RingType() != SugarToken::UNKNOWN_RING &&
      pdb_glycam->second.RingType() != sugarInfo.RingType())
  {
    mprintf("Warning: '%s' detected ring type is %s but ring type based on name is %s.\n",
            sugarName.c_str(), sugarInfo.ringStr(),
            SugarToken::ringStr(pdb_glycam->second.RingType()) );
    if (sugarInfo.RingType() != SugarToken::UNKNOWN_RING)
      resStatIn[rnum] = ResStatArray::SUGAR_NAME_MISMATCH;
    if (useSugarName_) {
      mprintf("Setting ring type bsaed on residue name.\n");
      sugarInfo.SetRingType( pdb_glycam->second.RingType() );
    } else if (sugar.IsMissingAtoms()) {
      mprintf("Warning: Residue was missing atoms; setting ring type based on residue name.\n");
      sugarInfo.SetRingType( pdb_glycam->second.RingType() );
    }
  }
    
  // Get glycam form string
  std::string formStr;
  if (sugarInfo.RingType() == SugarToken::PYRANOSE) {
    if (sugarInfo.Form() == SugarToken::ALPHA)
      formStr = "A";
    else
      formStr = "B";
  } else if (sugarInfo.RingType() == SugarToken::FURANOSE) {
    if (sugarInfo.Form() == SugarToken::ALPHA)
      formStr = "D";
    else
      formStr = "U";
  }
  // Sanity check
  if (formStr.empty()) {
    mprinterr("Internal Error: Could not set anomer type string.\n");
    return 1;
  }

  mprintf("\t  %s detected anomer type is %s(%s)-%s-%s\n", sugarName.c_str(),
         sugarInfo.formStr(), formStr.c_str(), sugarInfo.chirStr(),
         sugarInfo.ringStr());

  // Identify linkages to other residues.
  std::string linkcode = DetermineSugarLinkages(sugar, cmask, topIn, resStatIn,
                                                sugarBondsToRemove);
  if (linkcode.empty()) {
    //resStat_[rnum] = SUGAR_UNRECOGNIZED_LINKAGE;
    mprintf("Warning: Determination of sugar linkages failed.\n");
    return 0;
  }

  // Change PDB names to Glycam ones
  std::string resCode = pdb_glycam->second.GlycamCode();
  if (!hasGlycam_) {
    if (ChangePdbAtomNamesToGlycam(resCode, res, topIn, sugarInfo.Form()))
    {
      mprinterr("Error: Changing PDB atom names to Glycam failed.\n");
      return 1;
    }
  }

  // Modify residue char to indicate L form if necessary.
  // We do this here and not above so as not to mess with the
  // linkage determination.
  if (sugarInfo.Chirality() == SugarToken::IS_L) {
    for (std::string::iterator it = resCode.begin(); it != resCode.end(); ++it)
      *it = tolower( *it );
  }
  // Set new residue name
  NameType newResName( linkcode + resCode + formStr );
  if (!hasGlycam_) {
    mprintf("\t  Changing %s to Glycam resname: %s\n",
      topIn.TruncResNameOnumId(rnum).c_str(), *newResName);
    ChangeResName(res, newResName);
  } else {
    if (newResName.Truncated() != resName)
      mprintf("Warning: Detected glycam name '%s' differs from original name '%s'\n",
              *newResName, resName.c_str());
  }
  if (resStatIn[rnum] == ResStatArray::UNKNOWN)
    resStatIn[rnum] = ResStatArray::VALIDATED;
  return 0;
}

// -----------------------------------------------------------------------------
/** Attempt to find any missing linkages to the anomeric carbon in sugar. */
int SugarBuilder::FindSugarC1Linkages(int rnum1, int c_beg,
                                      Topology& topIn, Frame const& frameIn,
                                      NameType const& solventResName)
const
{
  Residue const& res1 = topIn.SetRes(rnum1);
  // If the anomeric atom is already bonded to another residue, skip this.
  for (Atom::bond_iterator bat = topIn[c_beg].bondbegin();
                           bat != topIn[c_beg].bondend(); ++bat)
  {
    if (topIn[*bat].ResNum() != rnum1) {
      if (debug_ > 0)
        mprintf("\tSugar %s anomeric carbon is already bonded to another residue, skipping.\n",
                topIn.TruncResNameOnumId(rnum1).c_str());
      return 0;
    }
  }
  // TODO pass these in
  // residue first atom to residue first atom cutoff^2
  const double rescut2 = 64.0;
  // bond cutoff offset
  const double offset = 0.2;
  // index of atom to be bonded to c_beg
  int closest_at = -1;
  // distance^2 of atom to be bonded to c_beg
  double closest_d2 = -1.0;

  Atom::AtomicElementType a1Elt = topIn[c_beg].Element(); // Should always be C
  if (debug_ > 0)
    mprintf("DEBUG: Anomeric ring carbon: %s\n", topIn.ResNameNumAtomNameNum(c_beg).c_str());
  // Loop over other residues
  for (int rnum2 = 0; rnum2 < topIn.Nres(); rnum2++)
  {
    if (rnum2 != rnum1) {
      Residue const& res2 = topIn.Res(rnum2);
      // Ignore solvent residues
      if (res2.Name() != solventResName) {
        int at1 = res1.FirstAtom();
        int at2 = res2.FirstAtom();
        // Initial residue-residue distance based on first atoms in each residue
        double dist2_1 = DIST2_NoImage( frameIn.XYZ(at1), frameIn.XYZ(at2) );
        if (dist2_1 < rescut2) {
          if (debug_ > 1)
            mprintf("DEBUG: Residue %s to %s = %f\n",
                    topIn.TruncResNameOnumId(rnum1).c_str(), topIn.TruncResNameOnumId(rnum2).c_str(),
                    sqrt(dist2_1));
          // Do the rest of the atoms in res2 to the anomeric carbon
          for (; at2 != res2.LastAtom(); ++at2)
          {
            if (!topIn[c_beg].IsBondedTo(at2)) {
              double D2 = DIST2_NoImage( frameIn.XYZ(c_beg), frameIn.XYZ(at2) );
              Atom::AtomicElementType a2Elt = topIn[at2].Element();
              double cutoff2 = Atom::GetBondLength(a1Elt, a2Elt) + offset;
              cutoff2 *= cutoff2;
              if (D2 < cutoff2) {
                if (debug_ > 1)
                  mprintf("DEBUG: Atom %s to %s = %f\n",
                          topIn.AtomMaskName(c_beg).c_str(), topIn.AtomMaskName(at2).c_str(), sqrt(D2));
                if (closest_at == -1) {
                  closest_at = at2;
                  closest_d2 = D2;
                } else if (D2 < closest_d2) {
                  mprintf("\t  Atom %s (%f Ang.) is closer than %s (%f Ang.).\n",
                          topIn.ResNameNumAtomNameNum(at2).c_str(), sqrt(D2),
                          topIn.ResNameNumAtomNameNum(closest_at).c_str(), sqrt(closest_d2));
                  closest_at = at2;
                  closest_d2 = D2;
                }
                //mprintf("\t  Adding bond between %s and %s\n",
                //        topIn.ResNameNumAtomNameNum(c_beg).c_str(),
                //        topIn.ResNameNumAtomNameNum(at2).c_str());
                //topIn.AddBond(c_beg, at2);
              }
            }
          } // END loop over res2 atoms
        } // END res1-res2 distance cutoff
      } // END res2 is not solvent
    } // END res1 != res2
  } // END res2 loop over other residues
  if (closest_at != -1) {
    mprintf("\t  Adding bond between %s and %s\n",
            topIn.ResNameNumAtomNameNum(c_beg).c_str(),
            topIn.ResNameNumAtomNameNum(closest_at).c_str());
    topIn.AddBond(c_beg, closest_at);
  }
  return 0;
}

/** Try to fix issues with sugar structure before trying to identify. */
int SugarBuilder::FixSugarsStructure(std::vector<Sugar>& sugarResidues,
                                            std::string const& sugarMaskStr,
                                            Topology& topIn, Frame& frameIn,
                                            bool c1bondsearch, bool splitres,
                                            NameType const& solventResName) 
const
{
  sugarResidues.clear();
  AtomMask sugarMask(sugarMaskStr);
  mprintf("\tLooking for sugars selected by '%s'\n", sugarMask.MaskString());
  if (topIn.SetupIntegerMask( sugarMask )) return 1;
  //sugarMask.MaskInfo();
  mprintf("\tSelected %i sugar atoms.\n", sugarMask.Nselected());
  if (sugarMask.None()) {
    mprintf("Warning: No sugar atoms selected by %s\n", sugarMask.MaskString());
    return 0;
  }
  //CharMask cmask( sugarMask.ConvertToCharMask(), sugarMask.Nselected() );
  // Get sugar residue numbers.
  Iarray sugarResNums = topIn.ResnumsSelectedBy( sugarMask );
  // Try to identify sugar rings.
  for (Iarray::const_iterator rnum = sugarResNums.begin();
                              rnum != sugarResNums.end(); ++rnum)
  {
    Sugar sugar = IdSugarRing(*rnum, topIn);
    if (sugar.Status() != Sugar::SETUP_OK) {
      mprintf("Warning: Problem identifying atoms for sugar '%s'\n",
              topIn.TruncResNameOnumId(*rnum).c_str());
    }
/*
    if (stat == ID_ERR) {
      if (errorsAreFatal_) {
        mprinterr("Error: Problem identifying sugar ring for %s\n",
                  topIn.TruncResNameOnumId(*rnum).c_str());
        return 1;
      } else
        mprintf("Warning: Problem identifying sugar ring for %s\n",
                topIn.TruncResNameOnumId(*rnum).c_str());
    }*/
    //if (!sugar.NotSet()) {
      sugarResidues.push_back( sugar );
      if (debug_ > 0)
        sugarResidues.back().PrintInfo(topIn);
    //}
  }

  if (c1bondsearch) {
    // Loop over sugar indices to see if anomeric C is missing bonds
    for (std::vector<Sugar>::const_iterator sugar = sugarResidues.begin();
                                            sugar != sugarResidues.end(); ++sugar)
    {
      if (sugar->NotSet()) continue;
      int anomericAtom = sugar->AnomericAtom();
      int rnum = sugar->ResNum(topIn);
      if (FindSugarC1Linkages(rnum, anomericAtom, topIn, frameIn, solventResName)) {
        mprinterr("Error: Search for bonds to anomeric carbon '%s' failed.\n",
                  topIn.AtomMaskName(anomericAtom).c_str());
        return 1;
      }
    }
  }
  //DEBUG
  //for (std::vector<Sugar>::const_iterator sugar = sugarResidues.begin();
  //                                    sugar != sugarResidues.end(); ++sugar)
  //  sugar->PrintInfo(topIn);
  FxnGroupBuilder FGB(debug_);
  if (splitres) {
    // Loop over sugar indices to see if residues have ROH that must be split off
    for (std::vector<Sugar>::iterator sugar = sugarResidues.begin();
                                      sugar != sugarResidues.end(); ++sugar)
    {
      if (sugar->NotSet()) continue;
      if (FGB.CheckIfSugarIsTerminal(*sugar, topIn, frameIn)) {
        mprinterr("Error: Checking if sugar %s has terminal functional groups failed.\n",
                  topIn.TruncResNameOnumId(sugar->ResNum(topIn)).c_str());
        return 1;
      }
    } // End loop over sugar indices

    // Loop over chain indices to see if residues need to be split
    for (std::vector<Sugar>::iterator sugar = sugarResidues.begin();
                                      sugar != sugarResidues.end(); ++sugar)
    {
      if (sugar->NotSet()) continue;
      if (FGB.CheckForFunctionalGroups(*sugar, topIn, frameIn)) {
        mprinterr("Error: Checking if sugar %s has functional groups failed.\n",
                 topIn.TruncResNameOnumId( sugar->ResNum(topIn) ).c_str());
        return 1;
      }
    }
    //DEBUG
    //for (std::vector<Sugar>::const_iterator sugar = sugarResidues.begin();
    //                                    sugar != sugarResidues.end(); ++sugar)
    //  sugar->PrintInfo(topIn);
  }


  return 0;
}

