#include <cstdio>
#include "../../../src/ArgList.h"

/** Convert CSV containing Abbreviation,Name,PDB,Glycam to
  * DAT file containing #Name GLYCAM PDB
  */
int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr,"Provide CSV file\n");
    return 1;
  }

  FILE* infile = fopen( argv[1], "rb" );
  if (infile == 0) return 1;

  char buffer[1024];
  static const int num = 1023;

  char* line = 0;
  int lineNum = 1;
  int err = 0;
  while ( (line = fgets( buffer, num, infile )) != 0 )
  {
    //printf("%s", line);
    ArgList argline(line, ",");
    if (argline.Nargs() != 4) {
      fprintf(stderr,"Error: Line %i has %i columns separated by commas, expected 4.\n", argline.Nargs());
      fprintf(stderr,"%s", line);
      err = 1;
      break;
    }
    lineNum++;
  }
  fclose( infile );

  return err;
}
