#include "SugarBuilder.h"
#include "../CpptrajStdio.h"

using namespace Cpptraj::Structure;

/** Load reduced internal PDB to Glycam map. */
void SugarBuilder::SetGlycamPdbResMap() {
  pdb_to_glycam_.insert( PairType("64K",
    SugarToken("alpha-D-arabinopyranose", "A", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("AHR",
    SugarToken("alpha-L-arabinofuranose", "A", ALPHA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("ARA",
    SugarToken("alpha-L-arabinopyranose", "A", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("ARB",
    SugarToken("beta-L-arabinopyranose", "A", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("AXR",
    SugarToken("methyl alpha-D-arabinofuranoside", "A", ALPHA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("BXX",
    SugarToken("beta-D-arabinofuranose", "A", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("BXY",
    SugarToken("alpha-D-arabinofuranose", "A", ALPHA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("FUB",
    SugarToken("beta-L-arabinofuranose", "A", BETA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("SEJ",
    SugarToken("beta-D-arabinopyranose", "A", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("LDY",
    SugarToken("alpha-D-lyxopyranose", "D", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z4W",
    SugarToken("beta-D-lyxopyranose", "D", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("0MK",
    SugarToken("beta-L-ribopyranose", "R", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("32O",
    SugarToken("beta-L-ribofuranose", "R", BETA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("BDR",
    SugarToken("beta-D-ribofuranose", "R", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("RIB",
    SugarToken("alpha-D-ribofuranose", "R", ALPHA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("RIP",
    SugarToken("beta-D-ribopyranose", "R", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("YYM",
    SugarToken("alpha-D-ribopyranose", "R", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z6J",
    SugarToken("alpha-L-ribofuranose", "R", ALPHA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("HSY",
    SugarToken("alpha-L-xylopyranose", "X", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("HSZ",
    SugarToken("beta-D-xylopyranose", "X", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("LXC",
    SugarToken("beta-L-xylopyranose", "X", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("XYP",
    SugarToken("beta-D-xylopyranose", "X", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("XYS",
    SugarToken("alpha-D-xylopyranose", "X", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("XYZ",
    SugarToken("beta-D-xylofuranose", "X", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("AFD",
    SugarToken("alpha-D-allopyranose", "N", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("ALL",
    SugarToken("beta-D-allopyranose", "N", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("VDS",
    SugarToken("beta-D-allofuranose", "N", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("VDV",
    SugarToken("alpha-D-allofuranose", "N", ALPHA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("WOO",
    SugarToken("beta-L-allopyranose", "N", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z2D",
    SugarToken("alpha-L-allopyranose", "N", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("3MK",
    SugarToken("beta-L-altropyranose", "E", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("SHD",
    SugarToken("alpha-D-altropyranose", "E", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z6H",
    SugarToken("alpha-L-altropyranose", "E", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GAL",
    SugarToken("beta-D-galactopyranose", "L", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GIV",
    SugarToken("beta-L-galactopyranose", "L", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GLA",
    SugarToken("alpha-D-galactopyranose", "L", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GXL",
    SugarToken("alpha-L-galactopyranose", "L", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GZL",
    SugarToken("beta-D-galactofuranose", "L", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("BGC",
    SugarToken("beta-D-glucopyranose", "G", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GLC",
    SugarToken("alpha-D-glucopyranose", "G", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GU4",
    SugarToken("2,3,4,6-tetra-O-sulfonato-alpha-D-glucopyranose", "G", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("MGL",
    SugarToken("methyl beta-D-glucopyranoside", "G", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z8T",
    SugarToken("beta-L-glucopyranose", "G", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("4GL",
    SugarToken("alpha-D-gulopyranose", "K", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GL0",
    SugarToken("beta-D-gulopyranose", "K", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GUP",
    SugarToken("alpha-L-gulopyranose", "K", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z8H",
    SugarToken("beta-L-gulopyranose", "K", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("BMA",
    SugarToken("beta-D-mannopyranose", "M", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("MAN",
    SugarToken("alpha-D-mannopyranose", "M", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("A5C",
    SugarToken("alpha-L-talofuranose", "T", ALPHA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("SDY",
    SugarToken("beta-D-talopyranose", "T", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("ZEE",
    SugarToken("beta-L-talopyranose", "T", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("BDF",
    SugarToken("beta-D-fructopyranose", "C", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("FRU",
    SugarToken("beta-D-fructofuranose", "C", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("LFR",
    SugarToken("beta-L-fructofuranose", "C", BETA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("YYJ",
    SugarToken("1,3,4,6-tetra-O-sulfo-beta-D-fructofuranose", "C", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("Z9N",
    SugarToken("alpha-D-fructofuranose", "C", ALPHA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("PSV",
    SugarToken("alpha-D-psicofuranose", "P", ALPHA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("SF6",
    SugarToken("alpha-L-psicofuranose", "P", ALPHA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("SF9",
    SugarToken("beta-L-psicofuranose", "P", BETA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("TTV",
    SugarToken("beta-D-psicofuranose", "P", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("SOE",
    SugarToken("alpha-L-sorbopyranose", "B", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("UEA",
    SugarToken("beta-D-sorbofuranose", "B", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("T6T",
    SugarToken("alpha-D-tagatopyranose", "J", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("FCA",
    SugarToken("alpha-D-fucopyranose", "F", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("FCB",
    SugarToken("beta-D-fucopyranose", "F", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("FUC",
    SugarToken("alpha-L-fucopyranose", "F", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("FUL",
    SugarToken("beta-L-fucopyranose", "F", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GYE",
    SugarToken("beta-D-fucofuranose", "F", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("MXY",
    SugarToken("2-O-methyl-beta-L-fucopyranose", "F", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("MXZ",
    SugarToken("2-O-methyl-alpha-L-fucopyranose", "F", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("G6D",
    SugarToken("alpha-D-quinovopyranose", "Q", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("YYK",
    SugarToken("beta-D-quinovopyranose", "Q", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("RAM",
    SugarToken("alpha-L-rhamnopyranose", "H", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("RM4",
    SugarToken("beta-L-rhamnopyranose", "H", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("XXR",
    SugarToken("alpha-D-rhamnopyranose", "H", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("ADA",
    SugarToken("alpha-D-galactopyranuronic acid", "O", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GTR",
    SugarToken("beta-D-galactopyranuronic acid", "O", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("BDP",
    SugarToken("beta-D-glucopyranuronic acid", "Z", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GCU",
    SugarToken("alpha-D-glucopyranuronic acid", "Z", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GCV",
    SugarToken("4-O-methyl-alpha-D-glucopyranuronic acid", "Z", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("IDR",
    SugarToken("alpha-L-idopyranuronic acid", "U", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("IDS",
    SugarToken("2-O-sulfo-alpha-L-idopyranuronic acid", "U", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("A2G",
    SugarToken("2-acetamido-2-deoxy-alpha-D-galactopyranose", "V", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("ASG",
    SugarToken("2-acetamido-2-deoxy-4-O-sulfo-beta-D-galactopyranose", "V", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("NG6",
    SugarToken("2-acetamido-2-deoxy-6-O-sulfo-beta-D-galactopyranose", "V", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("NGA",
    SugarToken("2-acetamido-2-deoxy-beta-D-galactopyranose", "V", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("YYQ",
    SugarToken("2-acetamido-2-deoxy-alpha-L-galactopyranose", "V", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("NAG",
    SugarToken("2-acetamido-2-deoxy-beta-D-glucopyranose", "Y", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("NDG",
    SugarToken("2-acetamido-2-deoxy-alpha-D-glucopyranose", "Y", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("NGZ",
    SugarToken("2-acetamido-2-deoxy-alpha-L-glucopyranose", "Y", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("BM3",
    SugarToken("2-acetamido-2-deoxy-alpha-D-mannopyranose", "W", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("BM7",
    SugarToken("2-acetamido-2-deoxy-beta-D-mannopyranose", "W", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("SIA",
    SugarToken("N-acetyl-alpha-neuraminic acid", "S", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("SLB",
    SugarToken("N-acetyl-beta-neuraminic acid", "S", BETA, IS_D, PYRANOSE)) );
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
int Exec_PrepareForLeap::LoadGlycamPdbResMap(std::string const& fnameIn)
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
