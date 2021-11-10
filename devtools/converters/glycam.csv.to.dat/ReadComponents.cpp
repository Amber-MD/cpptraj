#include "ArgList.h"
#include "CIFfile.h"
#include "GlycamPdbResMap.h"

/** Extract the "id" field from data with given name. */
std::string ExtractId(CIFfile const& cif, std::string const& name) {
  std::string id;
  printf("NAME: %s", name.c_str());
  CIFfile::DataBlock const& block = cif.GetBlockWithColValue("_chem_comp", "name", name);
  if (!block.empty()) {
    //block.ListData();
    id.assign( block.Data("id") );
    std::string replaced_by = block.Data("pdbx_replaced_by");
    if (!replaced_by.empty() && replaced_by != "?")
      printf(" 'replaced by %s' ", replaced_by.c_str());
  }
  return id;
}

/** Extract the "name" field from data with given id. */
std::string ExtractName(CIFfile const& cif, std::string const& id) {
  std::string name;
  //printf("ID: %s", id.c_str());
  CIFfile::DataBlock const& block = cif.GetBlockWithColValue("_chem_comp", "id", id);
  if (!block.empty()) {
    name.assign( block.Data("name") );
  } else {
    // Try pdbx_replaced_by (XYP only?)
    CIFfile::DataBlock const& block2 = cif.GetBlockWithColValue("_chem_comp", "pdbx_replaced_by", id);
    if (!block2.empty())
      name.assign( block2.Data("name") );
  }
  return name;
}

void NameToForm(std::string const& name, char glyc,
                std::string& form, std::string& chirality, std::string& ring)
{
  // alpha-L-arabinopyranose
  ArgList list1(name, "- ");
  if (list1.hasKey("alpha"))
    form = "A";
  else if (list1.hasKey("beta"))
    form = "B";
  else
    fprintf(stdout,"ERROR: Could not extract form.\n");
  if (glyc == 'S')
    // S is always D for glycam
    chirality = "D";
  else if (list1.hasKey("D"))
    chirality = "D";
  else if (list1.hasKey("L"))
    chirality = "L";
  else
    fprintf(stdout,"ERROR: Could not extract chirality.\n");
  ArgList remain = list1.RemainingArgs();
  ring = "P";
  for (int iarg = 0; iarg != remain.Nargs(); iarg++)
  {
    std::size_t found = remain[iarg].find("furanose");
    if (found != std::string::npos) {
      ring = "F";
      break;
    }
    //else
    ////found = list1[2].find("furanose")
    //// assume P
    //ring = "P";
  }
}

/** Print ID to stdout. */
void PrintId(std::string const& name, std::string const& id, std::string const& glyc) {
  if (!id.empty()) {
    printf(" %s %s %s", name.c_str(), id.c_str(), glyc.c_str());
    //std::string form, chirality, ring;
    //NameToForm(name, glyc, form, chirality, ring);
    //fprintf(direct_pdb_to_glycam, "%s %s %s %s %s\n",
    //        id.c_str(), glyc.c_str(), form.c_str(), chirality.c_str(), ring.c_str());
  }
  printf("\n\n");
}

typedef std::vector<std::string> Sarray;

/** Loop over alpha, beta, D, and L forms. */
void Loops(CIFfile const& cif, std::string const& lead,
           Sarray const& prefixes, Sarray const& suffixes, Sarray const& glycam)
{
  Sarray forms;
  forms.push_back("alpha");
  forms.push_back("beta");

  Sarray chirality;
  chirality.push_back("D");
  chirality.push_back("L");

  Sarray::const_iterator glyc = glycam.begin();
  for (Sarray::const_iterator pref = prefixes.begin(); pref != prefixes.end(); ++pref, ++glyc) {
    for (Sarray::const_iterator suff = suffixes.begin(); suff != suffixes.end(); ++suff) {
      for (Sarray::const_iterator form = forms.begin(); form != forms.end(); ++form) {
        for (Sarray::const_iterator chir = chirality.begin(); chir != chirality.end(); ++chir) {
          std::string name;
          if (!lead.empty())
            name.assign( lead + "-" );
          name.append( *form + "-" + *chir + "-" + *pref + *suff );
          std::string id = ExtractId(cif, name);
          PrintId( name, id, *glyc );
        }
      }
    }
  }
}

/** Loops, no lead string. */
void Loops(CIFfile const& cif, Sarray const& prefixes, Sarray const& suffixes,
           Sarray const& glycam)
{
  Loops(cif, std::string(""), prefixes, suffixes, glycam);
}

/** Get PDB codes from sugar names */
int GetPdbCodesFromNames(CIFfile& cif) {
  Sarray prefixes; Sarray glycam;
  prefixes.push_back("arabino"); glycam.push_back("A"); // A
  prefixes.push_back("lyxo"); glycam.push_back("D"); // D
  prefixes.push_back("ribo"); glycam.push_back("R"); // R
  prefixes.push_back("xylo"); glycam.push_back("X"); // X
  prefixes.push_back("allo"); glycam.push_back("N"); // N
  prefixes.push_back("altro"); glycam.push_back("E"); // E
  prefixes.push_back("galacto"); glycam.push_back("L"); // L
  prefixes.push_back("gluco"); glycam.push_back("G"); // G
  prefixes.push_back("gulo"); glycam.push_back("K"); // K
  prefixes.push_back("ido"); glycam.push_back("I"); // I (glycam does not have?)
  prefixes.push_back("manno"); glycam.push_back("M"); // M
  prefixes.push_back("talo"); glycam.push_back("T"); // T
  prefixes.push_back("fructo"); glycam.push_back("C"); // C
  prefixes.push_back("psico"); glycam.push_back("P"); // P
  prefixes.push_back("sorbo"); glycam.push_back("B"); // B
  prefixes.push_back("tagato"); glycam.push_back("J"); // J
  prefixes.push_back("fuco"); glycam.push_back("F"); // F
  prefixes.push_back("quinovo"); glycam.push_back("Q"); // Q
  prefixes.push_back("rhamno"); glycam.push_back("H"); // H

  Sarray suffixes;
  suffixes.push_back("pyranose");
  suffixes.push_back("furanose");

  Loops(cif, prefixes, suffixes, glycam);

  // -----------------------------------
  Sarray pprefixes; Sarray pglycam;
  pprefixes.push_back("galacto"); pglycam.push_back("O"); // O
  pprefixes.push_back("gluco"); pglycam.push_back("Z"); // Z
  pprefixes.push_back("ido"); pglycam.push_back("U"); // U

  Sarray psuffixes;
  psuffixes.push_back("pyranuronic acid");

  Loops(cif, pprefixes, psuffixes, pglycam);

  // -----------------------------------
  std::string lead("2-acetamido-2-deoxy");
  Sarray prefixes3; Sarray glycam3;
  prefixes3.push_back("galacto"); glycam3.push_back("V"); // V
  prefixes3.push_back("gluco"); glycam3.push_back("Y"); // Y
  prefixes3.push_back("manno"); glycam3.push_back("W"); // W

  Sarray suffixes3;
  suffixes3.push_back("pyranose");

  Loops(cif, lead, prefixes3, suffixes3, glycam3);
  // -----------------------------------
  PrintId("N-acetyl-alpha-neuraminic acid", ExtractId(cif, "N-acetyl-alpha-neuraminic acid"), "S"); // SA
  PrintId("N-acetyl-beta-neuraminic acid", ExtractId(cif, "N-acetyl-beta-neuraminic acid"), "S"); // SB
  return 0;
}

/** Get PDB names from current Carbohydrate_PDB_Glycam_Names.txt file,
  * extract form/chirality.
  */
int CreateDirectPdbToGlycam(FILE* direct_pdb_to_glycam, CIFfile const& cif) {
  GlycamPdbResMap pdb_to_glycam;
  if (pdb_to_glycam.Load("../../../dat/Carbohydrate_PDB_Glycam_Names.txt")) return 1;

  printf("\t%zu entries in PDB to glycam name map.\n", pdb_to_glycam.size());
  printf("\tResidue name map:\n");
  for (GlycamPdbResMap::const_iterator mit = pdb_to_glycam.begin();
                                       mit != pdb_to_glycam.end(); ++mit)
  {
    printf("\t  %4s -> %c ", *(mit->first), mit->second);
    std::string name = ExtractName(cif, mit->first.Truncated());
    if (name.empty())
      name.assign("EMPTY");
    printf(" %s\n", name.c_str());
    std::string form, chirality, ring;
    NameToForm(name, mit->second, form, chirality, ring);
    fprintf(direct_pdb_to_glycam, "%s %c %s %s %s \"%s\"\n",
            *(mit->first), mit->second, form.c_str(), chirality.c_str(), ring.c_str(), name.c_str());

  }

  return 0;
}

/** Read everything from the given CIF file. */
int ReadCIF(const char* fname) {


  CIFfile cif;
  if (cif.Read(fname, 0)) return 1;


  // DEBUG
//  //CIFfile::DataBlock lastBlock = cif.GetDataBlock("_chem_comp");
//  CIFfile::DataBlock const& lastBlock = cif.GetBlockWithColValue("_chem_comp",
//                                                                "name",
//                                                                "alpha-D-arabinopyranose");
//  lastBlock.ListData();

  //if (GetPdbCodesFromNames(cif)) return 1;

  FILE* direct_pdb_to_glycam = fopen("direct_pdb_to_glycam.dat", "wb");
  if (direct_pdb_to_glycam == 0) {
    fprintf(stderr,"Error: Could not open direct pdb to glycam file.\n");
    return 1;
  }
  int err = CreateDirectPdbToGlycam(direct_pdb_to_glycam, cif);
  fclose(direct_pdb_to_glycam);
  return err;
}

int main(int argc, char** argv) {
  if (ReadCIF("components.cif")) return 1;
  return 0;
}
