#include <cstdio>
#include <string>
#include <set>
#include <map>

using namespace std;

/// List of strings
typedef set<string> Slist;

/// Pair string to list of strings
typedef pair<string, Slist> Spair;

/// Map string to list of strings
typedef map<string, Slist> Smap;

/// Hold dependencies for source files
Smap Sources;

/// Hold dependencies for headers
Smap Headers;

enum FileType { SOURCE = 0, HEADER };

/** Add list of dependencies for the given file to map. */
void GetDependencies(string const& filename) {
  // Determine type
  FileType type;
  string ext;
  size_t found = filename.find_last_of(".");
  if (found != std::string::npos)
    ext = filename.substr(found);

  printf("FILE: %s  EXT: %s\n", filename.c_str(), ext.c_str());
}

// M A I N
int main(int argc, char** argv) {
  if (argc < 2) return  0;
  for (int i=1; i<argc; i++)
    GetDependencies(argv[i]);

  return 0;
}
