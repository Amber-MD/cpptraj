#include <cstdio>
#include <string>
#include "../../../src/ArgList.h"
#include "../../../src/StringRoutines.h"

/** Convert CSV containing Abbreviation,Name,PDB,Glycam to
  * DAT file containing #Name GLYCAM PDB
  */
int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr,"Provide CSV file\n");
    return 1;
  }

  FILE* infile = fopen( argv[1], "rb" );
  if (infile == 0) {
    fprintf(stderr,"Could not open input file.\n");
    return 1;
  }
  FILE* outfile = fopen( "../../../dat/Carbohydrate_PDB_Glycam_Names.txt", "wb");
  if (outfile == 0) {
    fprintf(stderr,"Could not open output file.\n");
    return 1;
  }
  fprintf(outfile,"# This file contains the mapping from common PDB names to Glycam residue codes.\n");
  //fprintf(outfile,"# Information largely obtained from http://glycam.org/docs/othertoolsservice/2016/06/09/3d-snfg-list-of-residue-names/#PDB\n");
  fprintf(outfile,"# Information obtained from mining the PDB chemical database (components.cif).\n");
  fprintf(outfile,"# Last updated %s\n", TimeString().c_str());

  char buffer[1024];
  static const int num = 1023;

  char* line = 0;
  int lineNum = 1;
  int err = 0;
  while ( (line = fgets( buffer, num, infile )) != 0 )
  {
    //printf("%s", line);
    ArgList argline(line, ",\r\n");
    if (argline.Nargs() != 4) {
      fprintf(stderr,"Error: Line %i has %i columns separated by commas, expected 4.\n", argline.Nargs());
      fprintf(stderr,"%s", line);
      err = 1;
      break;
    }
    // Header is: Abbreviation,Name,PDB,Glycam
    if (argline[0] == "Abbreviation") continue;
    //         0            1    2   3
    // infile: Abbreviation Name PDB Glycam
    if (argline[2] == ".") {
      printf("\"%s\" has no PDB name(s), skipping.\n", argline[1].c_str());
    } else if (argline[3] == ".") {
      printf("\"%s\" has no Gycam name, skipping.\n", argline[1].c_str());
    } else {
      // outfile: Name Glycam PDB(commas)
      std::string name_with_commas = "\"" + argline[1] + "\"";
      fprintf(outfile, "%-30s %3s ", name_with_commas.c_str(), argline[3].c_str());
      ArgList pdbnames(argline[2], " ");
      for (int n = 0; n < pdbnames.Nargs(); n++) {
        if (n > 0)
          fprintf(outfile,",");
        fprintf(outfile, "%s", pdbnames[n].c_str());
      }
      fprintf(outfile,"\n");
    }
    lineNum++;
  }

  // Add name map section manually
  fprintf(outfile, "\n# PDB to glycam atom name maps\n");
  fprintf(outfile, "V,W,Y C7,C2N  O7,O2N  C8,CME\n");
  fprintf(outfile, "S     C10,C5N O10,O5N C11,CME\n");

  // Add linkage res name map section manually
  fprintf(outfile, "\n# PDB to glycame linkage residue name maps\n");
  fprintf(outfile, "SER OLS\n");
  fprintf(outfile, "THR OLT\n");
  fprintf(outfile, "HYP OLP\n");
  fprintf(outfile, "ASN NLN\n");

  fclose( infile );

  return err;
}
