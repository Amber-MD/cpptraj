/*
 * FindDepend2
 * Dan Roe 2018
 * Given a list of source files with include directives, print a list of
 * dependencies. Will ignore headers that are not directly a part of
 * CPPTRAJ.
 */
#include <cstdio>
#include <cctype>
#include <cstring>
#include <string>
#include <set>
#include <map>
#include <vector>

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
  if (strcmp(headername,"mkl.h")==0) return true;
  if (strcmp(headername,"zlib.h")==0) return true;
  if (strcmp(headername,"bzlib.h")==0) return true;
  if (strcmp(headername,"netcdf.h")==0) return true;
  if (strcmp(headername,"pnetcdf.h")==0) return true;
  if (strcmp(headername,"sander.h")==0) return true;
  if (strcmp(headername,"omp.h")==0) return true;
  if (strcmp(headername,"OpenMM.h")==0) return true;
  if (strncmp(headername,"readline",8)==0) return true;
  if (strncmp(headername,"xdrfile",7)==0) return true;
  if (strncmp(headername,"lisp",4)==0) return true; // for readline
  return false;
}

/** Try to simplify a directory prefix that may contain '../'. */
std::string SimplifyDirPrefix(std::string const& dirPrefix) {
  if (dirPrefix.empty()) return dirPrefix;

  std::vector<long int> slashes;
  for (std::string::const_iterator it = dirPrefix.begin(); it != dirPrefix.end(); ++it)
    if ( *it == '/' ) {
      slashes.push_back( it - dirPrefix.begin() );
      //printf("\tSlash at %li\n", slashes.back());
    }
  // If only one slash, nothing to simplify
  if (slashes.size() < 2) return dirPrefix;
  //           012
  // Check for ../
  //           210
  for (std::vector<long int>::const_iterator it = slashes.begin(); it != slashes.end(); ++it)
  {
    if (*it > 1) {
      if ( dirPrefix[*it-2] == '.' && dirPrefix[*it-1] == '.' ) {
        //printf("\t../ at %li\n", *it);
        // If this is ../ corresponding to 2nd slash, negates everything from the beginning to here
        if (it - slashes.begin() == 1) {
          std::string newStr = dirPrefix.substr(*it+1, dirPrefix.size() - *it+1);
          //printf("\tNew prefix= %s\n", newStr.c_str());
          return newStr;
        }
      }
    }
  }
  return dirPrefix;
}

void SplitIntoDirAndBase(std::string& filename, std::string& dirPrefix, std::string& baseName) {
  // Determine path
  size_t found = filename.find_last_of("/");
  if (found == std::string::npos) {
    baseName = filename;
    //dirPrefix_.clear();
  } else {
    baseName = filename.substr(found+1);
    dirPrefix = filename.substr(0, found+1);
  }
  //printf("DEBUG: Dir='%s' Base='%s'\n", dirPrefix.c_str(), baseName.c_str());
  // Ascertain whether we can simplify the directory prefix
  std::string newDirPrefix = SimplifyDirPrefix( dirPrefix );
  if (!dirPrefix.empty() && newDirPrefix.empty()) {
    dirPrefix.clear();
    filename = baseName;
  }
}

/** Add list of dependencies for the given file to appropriate map. */
void GetDependencies(string const& filenameIn) {
  std::string filename = filenameIn;
  char buffer[BUFFERSIZE+1];
  char headername[BUFFERSIZE+1];
  // Determine path
  string baseName, dirPrefix;
  SplitIntoDirAndBase(filename, dirPrefix, baseName);

  // Determine type
  FileType type;
  string ext;
  size_t found = filename.find_last_of(".");
  if (found != string::npos)
    ext = filename.substr(found);

  //printf("FILE: %s  EXT: %s\n", filename.c_str(), ext.c_str());
  if (ext == ".cpp" || ext == ".c" || ext == ".cu") {
    type = SOURCE;
    // Each source file should only be accessed once
    Smap::iterator it = Sources.find( filename );
    if (it != Sources.end()) {
      fprintf(stderr,"Error: Source '%s' is being looked at more than once.\n", filename.c_str());
      return;
    }
  } else if (ext == ".h" || ext == ".cuh") {
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
          if (!IgnoreHeader(headername)) {
            std::string newFileName = dirPrefix + string(headername);
            std::string newdir, newbase;
            SplitIntoDirAndBase(newFileName, newdir, newbase);
            depends.insert( newFileName );
          }
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
