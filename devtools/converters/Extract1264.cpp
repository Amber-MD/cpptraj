#include "ArgList.h"
#include "BufferedLine.h"
#include <cstdio>

int main(int argc, char** argv)
{
  BufferedLine infile;
  if (infile.OpenFileRead(argv[1])) return 1;

  const char* ptr = infile.Line();
  while (ptr != 0) {
    ArgList line(ptr, "{}':, ");
    //line.PrintDebug();
    for (int iarg = 0; iarg < line.Nargs(); iarg += 2)
      printf("  c4params_.insert( NameMapPair(\"%s\", %s) );\n", line[iarg].c_str(), line[iarg+1].c_str());
    ptr = infile.Line();
  }

  infile.CloseFile();
  return 0;
}
