#include <cstdlib> // atoi
#include "TinkerFile.h"
#include "ArgList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
TinkerFile::TinkerFile() : natom_(0), hasBox_(false) {}

static inline int SetNatomAndTitle(ArgList& firstLine, int& natom, std::string& title) {
  if (firstLine.Nargs() != 2) return 1;
  natom = firstLine.getNextInteger( -1 );
  if (natom < 1) return 1;
  title = firstLine.GetStringNext();
  if (title.empty()) return 1;
  return 0;
}

bool TinkerFile::ID_Tinker(CpptrajFile& fileIn) {
  // NOTE: ASSUME FILE SET UP FOR READ
  if (fileIn.OpenFile()) return false;
  ArgList firstLine( fileIn.NextLine() );
  // First line should have <natom> <title> only
  int natom = 0;
  std::string title;
  bool isTinker = ( SetNatomAndTitle(firstLine, natom, title)==0 );
  fileIn.CloseFile();
  return isTinker;
}

int TinkerFile::OpenTinker() {
  if (tinkerName_.empty()) {
    mprinterr("Internal Error: Tinker file name not set.\n");
    return 1;
  }
  if (file_.OpenFileRead( tinkerName_ )) return 1;
  ArgList firstLine( file_.Line() );
  if ( SetNatomAndTitle(firstLine, natom_, title_) ) {
    mprinterr("Error: Could not get # atoms / title from Tinker file.\n");
    return 1;
  }
  // Are box coords present? If not, next line read in should be 1 followed by
  // atom info and the line after that should be 2.
  hasBox_ = false;
  const char* secondptr = file_.Line();
  if (secondptr == 0) {
    mprinterr("Error: Could not get first atom line of Tinker file.\n");
    return 1;
  }
  // Third line could be zero if only 1 atom and no box.
  const char* thirdptr = file_.Line();
  if (natom_ == 1) {
    // If a third line was read, check if it is another title line. If so,
    // no box coordinates.
    if (thirdptr != 0) {
      firstLine.SetList( std::string(thirdptr), " " );
      int natom2;
      std::string title2; // TODO: Check natom/title match?
      if (SetNatomAndTitle(firstLine, natom2, title2))
        hasBox_ = true;
    } // else no third line read, no box.
  } else {
    if (thirdptr == 0) {
      mprinterr("Error: Could not get second atom line of Tinker file.\n");
      return 1;
    }
    firstLine.SetList( std::string(thirdptr), " " );
    // If the third line contains atom 1 there are box coords.
    int atomIdx = firstLine.IntegerAt( 0 );
    if (atomIdx < 1) {
      mprinterr("Error: Third line contains invalid atom index.\n");
      mprinterr("Error: %s", thirdptr);
      return 1;
    }
    if (atomIdx == 1)
      hasBox_ = true;
  }
  // Close and reopen the file.
  file_.CloseFile();
  return file_.OpenFileRead( tinkerName_ );
}

int TinkerFile::CheckTitleLine() {
  file_.TokenizeLine(" ");
  int lineNatom = atoi( file_.NextToken() );
  if (lineNatom != natom_) {
    mprinterr("Error: Number of atoms in Tinker file changes from %i to %i\n",
              "Error: at line %i\n", natom_, lineNatom, file_.LineNumber());
    return 1;
  }
  return 0;
}

/** \return 0 if no more frames to read.
  * \return -1 if an error occurs.
  * \return 1 if more frames to read.
  */
int TinkerFile::NextTinkerFrame() {
  // Title line
  if (file_.Line() == 0) return 0;
  if (CheckTitleLine()) return -1;
  // Box line if necessary
  if (hasBox_) {
    if (file_.Line() == 0) return -1;
  }
  for (int atidx = 0; atidx < natom_; atidx++)
    if (file_.Line() == 0) return -1;
  return 1;
}

int TinkerFile::ReadNextTinkerFrame(double* Xptr, double* box) {
  // Title line
  if (file_.Line() == 0) return 0;
  if (CheckTitleLine()) return -1;
  // Box line
  if (hasBox_) {
    if (file_.Line() == 0) return -1;
    int nbox = file_.TokenizeLine(" ");
    if (nbox != 6) {
      mprinterr("Error: In Tinker file line %i expected 6 box coords, got %i\n", 
                file_.LineNumber(), nbox);
      return -1;
    }
    for (int b = 0; b != nbox; b++)
      box[b] = atof( file_.NextToken() );
  }
  // Coords
  for (int atidx = 0; atidx < natom_; atidx++) {
    if (file_.Line() == 0) return -1;
    int ncol = file_.TokenizeLine(" ");
    if (ncol < 5) {
      mprinterr("Error: In Tinker file line %i expected at least 5 columns for atom, got %i\n",
                file_.LineNumber(), ncol);
      return -1;
    }
    file_.NextToken(); // Atom index
    file_.NextToken(); // Atom name
    Xptr[0] = atof( file_.NextToken() ); // X
    Xptr[1] = atof( file_.NextToken() ); // Y
    Xptr[2] = atof( file_.NextToken() ); // Z
    Xptr += 3;
  }
  return 1;
}
     
