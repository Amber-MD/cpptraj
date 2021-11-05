#include "ArgList.h"
#include "CIFfile.h"

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

/** Print ID to stdout. */
void PrintId(std::string const& name, std::string const& id) {
  if (!id.empty()) {
    printf(" %s %s", name.c_str(), id.c_str());
  }
  printf("\n\n");
}

typedef std::vector<std::string> Sarray;

/** Loop over alpha, beta, D, and L forms. */
void Loops(CIFfile const& cif, std::string const& lead,
           Sarray const& prefixes, Sarray const& suffixes)
{
  Sarray forms;
  forms.push_back("alpha");
  forms.push_back("beta");

  Sarray chirality;
  chirality.push_back("D");
  chirality.push_back("L");

  for (Sarray::const_iterator pref = prefixes.begin(); pref != prefixes.end(); ++pref) {
    for (Sarray::const_iterator suff = suffixes.begin(); suff != suffixes.end(); ++suff) {
      for (Sarray::const_iterator form = forms.begin(); form != forms.end(); ++form) {
        for (Sarray::const_iterator chir = chirality.begin(); chir != chirality.end(); ++chir) {
          std::string name;
          if (!lead.empty())
            name.assign( lead + "-" );
          name.append( *form + "-" + *chir + "-" + *pref + *suff );
          std::string id = ExtractId(cif, name);
          PrintId( name, id );
        }
      }
    }
  }
}

/** Loops, no lead string. */
void Loops(CIFfile const& cif, Sarray const& prefixes, Sarray const& suffixes) {
  Loops(cif, std::string(""), prefixes, suffixes);
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

  Sarray prefixes;
  prefixes.push_back("arabino"); // A
  prefixes.push_back("lyxo"); // D
  prefixes.push_back("ribo"); // R
  prefixes.push_back("xylo"); // X
  prefixes.push_back("allo"); // N
  prefixes.push_back("altro"); // E
  prefixes.push_back("galacto"); // L
  prefixes.push_back("gluco"); // G
  prefixes.push_back("gulo"); // K
  prefixes.push_back("ido"); // I (glycam does not have?)
  prefixes.push_back("manno"); // M
  prefixes.push_back("talo"); // T
  prefixes.push_back("fructo"); // C
  prefixes.push_back("psico"); // P
  prefixes.push_back("sorbo"); // B
  prefixes.push_back("tagato"); // J
  prefixes.push_back("fuco"); // F
  prefixes.push_back("quinovo"); // Q
  prefixes.push_back("rhamno"); // H

  Sarray suffixes;
  suffixes.push_back("pyranose");
  suffixes.push_back("furanose");

  Loops(cif, prefixes, suffixes);

  // -----------------------------------
  Sarray pprefixes;
  pprefixes.push_back("galacto"); // O
  pprefixes.push_back("gluco"); // Z
  pprefixes.push_back("ido"); // U

  Sarray psuffixes;
  psuffixes.push_back("pyranuronic acid");

  Loops(cif, pprefixes, psuffixes);

  // -----------------------------------
  std::string lead("2-acetamido-2-deoxy");
  Sarray prefixes3;
  prefixes3.push_back("galacto"); // V
  prefixes3.push_back("gluco"); // Y
  prefixes3.push_back("manno"); // W

  Sarray suffixes3;
  suffixes3.push_back("pyranose");

  Loops(cif, lead, prefixes3, suffixes3);
  // -----------------------------------
  PrintId("N-acetyl-alpha-neuraminic acid", ExtractId(cif, "N-acetyl-alpha-neuraminic acid")); // SA
  PrintId("N-acetyl-beta-neuraminic acid", ExtractId(cif, "N-acetyl-beta-neuraminic acid")); // SB

  return 0;
}

int main(int argc, char** argv) {
  if (ReadCIF("components.cif")) return 1;
  return 0;
}
