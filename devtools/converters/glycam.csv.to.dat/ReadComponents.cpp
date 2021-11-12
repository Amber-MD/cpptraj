#include "ArgList.h"
#include "CIFfile.h"
#include "GlycamPdbResMap.h"
#include "StringRoutines.h"

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
    // Do furano to catch furanoside as well
    std::size_t found = remain[iarg].find("furano");
    if (found != std::string::npos) {
      ring = "F";
      break;
    }
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

/** Read residue names and glycam codes from file. Write the
  * Carbohydrate_PDB_Glycam_Names.txt file in the 'dat' directory.
  */
int GenerateFile(std::string const& f_glycamnames, std::string const& f_resnames,
                 std::string const& f_outfile, CIFfile const& cif)
{
  typedef std::pair<std::string, std::string> Spair;
  //typedef std::map<std::string, std::string> Smap;
  // Use a vector to keep the ordering in the GlycamNames.txt file.
  typedef std::vector<Spair> SpairArray;
  // First read sugar names and glycam codes from infile
  SpairArray GlycamCode_to_sugarName;
  BufferedLine glycamnames;
  if (glycamnames.OpenFileRead(f_glycamnames)) return 1;
  const char* ptr = glycamnames.Line();
  while (ptr != 0) {
    if (*ptr != '#') {
      ArgList line(ptr, " ");
      if (line.Nargs() != 2) {
        fprintf(stderr,"Error: Expected only 2 columns in glycam name file.\n");
        return 1;
      }
      //printf("DEBUG: %s %s\n", line[0].c_str(), line[1].c_str());
      GlycamCode_to_sugarName.push_back(Spair(line[1], line[0]));
    }
    ptr = glycamnames.Line();
  }
  glycamnames.CloseFile();
  for (SpairArray::const_iterator it = GlycamCode_to_sugarName.begin();
                                  it != GlycamCode_to_sugarName.end(); ++it)
    printf("%s -> %s\n", it->first.c_str(), it->second.c_str());

  // Now read sugar residue names and their glycam codes
  typedef std::vector<std::string> Sarray;
  typedef std::pair<std::string, Sarray> CodeSarrayPair;
  typedef std::map<std::string, Sarray> CodeSarrayMap;
  CodeSarrayMap GlycamCode_to_resNames;
  BufferedLine resnames;
  if (resnames.OpenFileRead(f_resnames)) return 1;
  ptr = resnames.Line();
  while (ptr != 0) {
    if (*ptr != '#') {
      ArgList line(ptr, " ");
      if (line.Nargs() != 2) {
        fprintf(stderr,"Error: Expected only 2 columns in res name file.\n");
        return 1;
      }
      // TODO ensure valid code?
      CodeSarrayMap::iterator it = GlycamCode_to_resNames.find( line[1] );
      if (it == GlycamCode_to_resNames.end()) {
        // New glycam code
        std::pair<CodeSarrayMap::iterator,bool> ret =
          GlycamCode_to_resNames.insert( CodeSarrayPair( line[1], Sarray() ) );
        it = ret.first;
      }
      it->second.push_back(line[0]);
    }
    ptr = resnames.Line();
  }

  // Output
  CpptrajFile outfile;
  if (outfile.OpenWrite(f_outfile)) return 1;
  outfile.Printf("# This file contains the mapping from common PDB names to Glycam residue codes.\n");
  outfile.Printf("# Information obtained from mining the PDB chemical database (components.cif).\n");
  outfile.Printf("# Last updated %s\n", TimeString().c_str());
  outfile.Printf("#ResName GlycamCode Form Chirality RingType \"Name\"\n");

  for (SpairArray::const_iterator kt = GlycamCode_to_sugarName.begin();
                                  kt != GlycamCode_to_sugarName.end(); ++kt)
  {
    CodeSarrayMap::const_iterator it = GlycamCode_to_resNames.find( kt->first );
    if (it == GlycamCode_to_resNames.end()) {
      fprintf(stderr,"Error: Glycam code '%s' not found in resNames.\n", kt->first.c_str());
      return 1;
    }
    std::string sname = "\"" + kt->second + "\"";
    printf("%-30s %3s ", sname.c_str(), it->first.c_str());
    for (Sarray::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
    {
      if (jt != it->second.begin()) printf(",");
      printf("%s", jt->c_str());
      // Find the name in components
      std::string fullname = ExtractName(cif, *jt);
      if (fullname.empty())
        fullname.assign("EMPTY");
      //printf(" %s\n", fullname.c_str());
      std::string form, chirality, ring;
      NameToForm(fullname, it->first[0], form, chirality, ring);
      outfile.Printf("%s %s %s %s %s \"%s\"\n",
                    jt->c_str(), it->first.c_str(),
                    form.c_str(), chirality.c_str(), ring.c_str(), fullname.c_str());
    }
    printf("\n");
  }

  // Add name map section manually
  outfile.Printf("\n# PDB to glycam atom name maps\n");
  outfile.Printf("V,W,Y C7,C2N  O7,O2N  C8,CME\n");
  outfile.Printf("S     C10,C5N O10,O5N C11,CME\n");
  outfile.Printf("H     C6,C6M,B\n");

  // Add linkage res name map section manually
  outfile.Printf("\n# PDB to glycam linkage residue name maps\n");
  outfile.Printf("SER OLS\n");
  outfile.Printf("THR OLT\n");
  outfile.Printf("HYP OLP\n");
  outfile.Printf("ASN NLN\n");

  outfile.CloseFile();
  return 0;
}
    

/** Read everything from the given CIF file. */
int ReadCIF(const char* fname) {
  int err = 0;

  CIFfile cif;
  if (cif.Read(fname, 0)) return 1;


  // DEBUG
//  //CIFfile::DataBlock lastBlock = cif.GetDataBlock("_chem_comp");
//  CIFfile::DataBlock const& lastBlock = cif.GetBlockWithColValue("_chem_comp",
//                                                                "name",
//                                                                "alpha-D-arabinopyranose");
//  lastBlock.ListData();

  //if (GetPdbCodesFromNames(cif)) return 1;

  if (GenerateFile("GlycamNames.txt", "../../../dat/ResNames.sugar.dat", "temp.dat", cif))
    return 1;

/*  FILE* direct_pdb_to_glycam = fopen("direct_pdb_to_glycam.dat", "wb");
  if (direct_pdb_to_glycam == 0) {
    fprintf(stderr,"Error: Could not open direct pdb to glycam file.\n");
    return 1;
  }
  int err = CreateDirectPdbToGlycam(direct_pdb_to_glycam, cif);
  fclose(direct_pdb_to_glycam);*/
  return err;
}

int main(int argc, char** argv) {
  if (ReadCIF("components.cif")) return 1;
  return 0;
}
