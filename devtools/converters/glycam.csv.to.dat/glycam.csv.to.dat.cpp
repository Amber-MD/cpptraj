#include <cstdio>

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
  while ( (line = fgets( buffer, num, infile )) != 0 )
  {
    printf("%s", line);
  }
  fclose( infile );

  return 0;
}
