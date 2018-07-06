#include <cstdio>
#include <cctype>
#include <cstring>
#include <string>
#include <set>
#include <map>

using namespace std;

#define BUFFERSIZE 1023

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

/** \return true if this header should be ignored. */
bool IgnoreHeader(const char* headername) {
  if (strcmp(headername,"mpi.h")==0) return true;
  if (strcmp(headername,"zlib.h")==0) return true;
  if (strcmp(headername,"bzlib.h")==0) return true;
  if (strcmp(headername,"netcdf.h")==0) return true;
  if (strcmp(headername,"pnetcdf.h")==0) return true;
  if (strcmp(headername,"sander.h")==0) return true;
  if (strcmp(headername,"omp.h")==0) return true;
  if (strncmp(headername,"readline",8)==0) return true;
  if (strncmp(headername,"xdrfile",7)==0) return true;
  return false;
}

/** Add list of dependencies for the given file to map. */
void GetDependencies(string const& filename) {
  char buffer[BUFFERSIZE+1];
  char headername[BUFFERSIZE+1];
  // Determine type
  FileType type;
  string ext;
  size_t found = filename.find_last_of(".");
  if (found != string::npos)
    ext = filename.substr(found);

  //printf("FILE: %s  EXT: %s\n", filename.c_str(), ext.c_str());
  if (ext == ".cpp" || ext == ".c") {
    type = SOURCE;
    // Each source file should only be accessed once
    Smap::iterator it = Sources.find( filename );
    if (it != Sources.end()) {
      fprintf(stderr,"Error: Source '%s' is being looked at more than once.\n", filename.c_str());
      return;
    }
  } else if (ext == ".h") {
    type = HEADER;
    // If this header was already looked at return
    Smap::iterator it = Headers.find( filename );
    if (it != Headers.end()) {
      //printf("\tSkipping already-seen header %s\n", filename.c_str());
      return;
    }
  } else // Ignore all others for now
    return;

  FILE* infile = fopen(filename.c_str(), "rb");
  if (infile == 0) {
    fprintf(stderr,"Error: Could not open '%s'\n", filename.c_str());
    return;
  }
  // Go through file and grab includes.
  Slist depends;
  while ( fgets(buffer, BUFFERSIZE, infile) != 0 )
  {
    char* ptr = buffer;
    // Skip leading whitespace.
    while ( isspace(*ptr) ) {
      ptr++;
      if (*ptr=='\n' || *ptr=='\0') break;
    }
    // First non-whitespace charcter must be '#'
    if (*ptr == '#') {
      // Check if this is an include line - if not, skip
      char* includePosition = strstr(ptr,"include");
      if ( includePosition != 0 ) {
        // Skip over system headers
        if ( strchr(ptr, '<') == 0 ) {
          // Get header name - assume it is second string - dont include first "
          sscanf(includePosition, "%*s \"%s", headername);
          // Get rid of last "
          size_t pos = strlen(headername);
          if (headername[pos-1]=='"') headername[pos-1]='\0';
          if (!IgnoreHeader(headername))
            depends.insert( string(headername) );
        } //else
          //printf("\tSkipping system header line: %s", buffer);
      }
    }
  }
  fclose( infile );
  //printf("  %s depends:", filename.c_str());
  //for (Slist::const_iterator it = depends.begin(); it != depends.end(); ++it)
  //  printf(" %s", it->c_str());
  //printf("\n");
  //pair<Smap::iterator, bool> ret;
  if (type == SOURCE)
    Sources.insert( Spair(filename, depends) );
  else
    Headers.insert( Spair(filename, depends) );
  for (Slist::const_iterator it = depends.begin(); it != depends.end(); ++it)
     GetDependencies( *it ); 
}

/** Recursive function for expanding dependencies. */
void FillDependencies(Slist& depends, Slist const& dlist) {
  for (Slist::const_iterator hdr = dlist.begin(); hdr != dlist.end(); ++hdr) {
    depends.insert( *hdr );
    Smap::iterator ret = Headers.find( *hdr );
    if (ret != Headers.end())
      FillDependencies(depends, ret->second);
  }
}

// M A I N
int main(int argc, char** argv) {
  if (argc < 2) return  0;
  // Read through all source files; get dependencies and header dependencies.
  for (int i=1; i<argc; i++)
    GetDependencies(argv[i]);
  // Write out source dependencies and fill in header dependencies.
  for (Smap::const_iterator src = Sources.begin(); src != Sources.end(); ++src)
  {
    Slist depends;
    FillDependencies(depends, src->second);
    size_t found = src->first.find_last_of(".");
    string objname = src->first.substr(0,found) + ".o";
    printf("%s : %s", objname.c_str(), src->first.c_str());
    for (Slist::const_iterator it = depends.begin(); it != depends.end(); ++it)
      printf(" %s", it->c_str());
    printf("\n");
  }

  return 0;
}
