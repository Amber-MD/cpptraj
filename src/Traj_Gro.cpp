#include <cstdio>
#include <cstdlib>
#include "Traj_Gro.h"
#include "StringRoutines.h"
#include "CpptrajStdio.h"

bool Traj_Gro::ID_TrajFormat(CpptrajFile& infile) {
  // Title line, atoms line, then resnum, resname, atomname, atomnum, X, Y, Z
  if (infile.OpenFile()) return false;
  if (infile.NextLine() == 0) return false; // Title
  if (infile.NextLine() == 0) return false; // Natom
  const char* ptr = infile.NextLine(); // First atom
  char resnum[6], resname[6], atname[6], atnum[6];
  float XYZ[3];
  int nread = sscanf(ptr, "%5c%5c%5c%5c%f %f %f", resnum, resname,
                     atname, atnum, XYZ, XYZ+1, XYZ+2);
  infile.CloseFile();
  return (nread == 7);
}

/** Assume time value cannot be negative. */
double Traj_Gro::GetTimeValue(const char* line) const {
  double timeVal = -1.0;
  if (line != 0 && *line != '\0') {
    const char* ptr = line;
    while ( *(ptr+2) != '\0' ) {
      if (*ptr == 't' && *(ptr+1) == '=') {
        timeVal = atof( ptr + 2 );
        break;
      }
      ++ptr;
    }
  }
  return timeVal;
}

int Traj_Gro::setupTrajin(std::string const& fname, Topology* trajParm)
{
  float fXYZ[9];
  // Open file for reading 
  if (file_.OpenFileRead( fname )) return TRAJIN_ERR;
  // Read the title. May contain time value, 't= <time>'
  const char* ptr = file_.Line();
  if (ptr == 0) {
    mprinterr("Error: Reading title.\n");
    return TRAJIN_ERR;
  }
  std::string title( ptr );
  bool hasTime = true;
  // TODO Is it OK to assume there will never be a negative time value?
  double timeVal = GetTimeValue( ptr );
  if (timeVal < 0.0) hasTime = false;
  /*std::string title = file_.GetLine();
  size_t pos = title.find("t=", 0);
  if (pos != std::string::npos) {
    std::string timeval = title.substr(pos + 2);
    if (!validDouble( timeval ))
      mprintf("Warning: Invalid time value: %s\n", timeval.c_str());
    else {
      time0_ = atof(timeval.c_str());
      hasTime = true;
    }
  } else
    time0_ = 0.0;*/
  // Read number of atoms
  ptr = file_.Line();
  if (ptr == 0) return TRAJIN_ERR;
  natom_ = atoi(ptr);
  if (natom_ < 1) {
    mprinterr("Error: Reading number of atoms.\n");
    return TRAJIN_ERR;
  }
  if (natom_ != trajParm->Natom()) {
    mprinterr("Error: Number of atoms %i does not match associated parm %s (%i)\n",
              natom_, trajParm->c_str(), trajParm->Natom());
    return TRAJIN_ERR;
  }
  // Read first atom to see if there are velocities
  int nread = sscanf(ptr, "%*5c%*5c%*5c%*5c%f %f %f %f %f %f",
                     fXYZ, fXYZ+1, fXYZ+2, fXYZ+3, fXYZ+4, fXYZ+5);
  bool hasV = false;
  if (nread == 6)
    hasV = true;
  else if (nread != 3) {
    mprinterr("Error: Reading first atom, expected 3 or 6 coordinates, got %i\n", nread);
    return TRAJIN_ERR;
  }
  // Read past the rest of the atoms
  for (int i = 1; i != natom_; i++)
    if (file_.Line() == 0) {
      mprinterr("Error: Reading atom %i of first frame.\n", i+1);
      return TRAJIN_ERR;
    }
    
  // Attempt to read box
  // v1(x) v2(y) v3(z) [v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)]
  ptr = file_.Line();
  Box groBox;
  if (ptr != 0) {
    nbox_ = sscanf( ptr, "%f %f %f %f %f %f %f %f %f", fXYZ, fXYZ+1, fXYZ+2,
                    fXYZ+3, fXYZ+4, fXYZ+5, fXYZ+6, fXYZ+7, fXYZ+8);
    if (nbox_ == 3) {
      // Box lengths. Convert nm->Ang. Assume orthogonal FIXME OK?
      double dXYZ[6];
      dXYZ[0] = (double)fXYZ[0];
      dXYZ[1] = (double)fXYZ[1];
      dXYZ[2] = (double)fXYZ[2];
      for (int i = 0; i != 3; i++) dXYZ[i] *= 10.0;
      dXYZ[3] = 90.0;
      dXYZ[4] = 90.0;
      dXYZ[5] = 90.0;
      groBox.SetBox( dXYZ );
    } else if (nbox_ == 9) {
      // Box vectors. Convert nm->Ang.
      Matrix_3x3 ucell;
      ucell[0] = (double)fXYZ[0]; ucell[1] = (double)fXYZ[3]; ucell[2] = (double)fXYZ[4]; // X
      ucell[3] = (double)fXYZ[5]; ucell[4] = (double)fXYZ[1]; ucell[5] = (double)fXYZ[6]; // Y
      ucell[6] = (double)fXYZ[7]; ucell[7] = (double)fXYZ[8]; ucell[8] = (double)fXYZ[2]; // Z
      for (int i = 0; i != 9; i++) ucell[i] *= 10.0;
      groBox.SetBox( ucell );
    } // Otherwise assume no box.
  }
  // Set trajectory information. No temperature info.
  SetCoordInfo( CoordinateInfo(groBox, hasV, false, hasTime) );
  SetTitle( title );
  // Check for multiple frames. If nothing was read above, 1 frame, no box. If
  // box info was read but nothing more can be read, 1 frame. Otherwise there
  // are more frames.
  bool hasMultipleFrames = false;
  if (ptr != 0) {
    if (groBox.Type() == Box::NOBOX) // Assume another title read.
      hasMultipleFrames = true;
    else {
      ptr = file_.Line(); // Read line after box info
      if (ptr != 0)
        hasMultipleFrames = true;
    }
  }
  int nframes = 1;
  if (hasMultipleFrames) {
    int linesToRead = natom_;
    if (groBox.Type() != Box::NOBOX)
      linesToRead += 1;
    // Since there is no guarantee that each frame is the same size we cannot
    // just seek. Blank reads for as many times as possible. Should currently
    // be positioned at the title line of the next frame.
    while (ptr != 0 ) {
      ptr = file_.Line(); // Natoms
      int Nat = atoi(ptr);
      if (Nat != natom_) {
        mprinterr("Error: Frame %i # atoms (%i) does not match first frame (%i).\n"
                  "Error: Only reading %i frames.\n", nframes+1, Nat, natom_, nframes);
        break;
      }
      for (int i = 0; i != linesToRead; i++)
        ptr = file_.Line();
      if (ptr == 0) break;
      nframes++;
      ptr = file_.Line(); // Next title or EOF
    }
  }
  file_.CloseFile();
  return nframes;
}
