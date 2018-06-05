#include <cstring> // strncmp
#include "CharmmParam.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"

int CharmmParam::ReadParams(std::string const& nameIn) {
  BufferedLine infile;

  mprintf("\tAttempting to read CHARMM parameters from '%s'\n", nameIn.c_str());

  if (infile.OpenFileRead( nameIn )) {
    mprinterr("Error: Could not open file '%s'\n", nameIn.c_str());
    return 1;
  }

  // NOTE: IGNORE is a debug option
  enum SectionType { UNKNOWN = 0, ATOMS, IGNORE };
  SectionType currentSection = UNKNOWN;

  const char* ptr = infile.Line();
  while (ptr != 0) {
    if (*ptr == '*') {
      // Title line.
      mprintf("DEBUG: Title: %s\n", ptr);
    } else if (*ptr == '!') {
      // Comment line.
      mprintf("DEBUG: Comment: %s\n", ptr);
    } else if (*ptr == '\0') {
      // A blank line denotes section end, put us back into the unknown
      mprintf("DEBUG: Resetting to UNKNOWN at line %i\n", infile.LineNumber());
      currentSection = UNKNOWN;
    } else {
      if (currentSection == UNKNOWN) {
        if (strncmp( ptr, "ATOMS", 5 ) == 0) {
          mprintf("DEBUG: Entering ATOMS section (%i).\n", infile.LineNumber());
          currentSection = ATOMS;
        } else {
          mprintf("DEBUG: Ignoring unknown section '%s' (%i)\n", ptr, infile.LineNumber());
          currentSection = IGNORE;
        }
      }
    }
    ptr = infile.Line();
  }

  return 0;
}
