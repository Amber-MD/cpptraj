#include <cstdio>
#include <cstdlib>
#include "Traj_Gro.h"
#include "Topology.h"
#include "Frame.h"
#include "CpptrajFile.h"
#include "StringRoutines.h"
#include "CpptrajStdio.h"
#include "Constants.h"

// Traj_Gro::ID_TrajFormat()
bool Traj_Gro::ID_TrajFormat(CpptrajFile& infile) {
  // Title line, atoms line, then resnum, resname, atomname, atomnum, X, Y, Z
  if (infile.OpenFile()) return false;
  int nread = 0;
  if (infile.NextLine() != 0) { // Title
    const char* ptr = infile.NextLine(); // Natom
    if (ptr != 0) {
      // Advance to first non-space character
      const char* aptr = ptr;
      while (*aptr == ' ') ++aptr;
      // Ensure only a single valid integer on # atoms line. This is done to
      // avoid confusion with Amber ASCII coords.
      std::string natom_str( aptr );
      RemoveTrailingWhitespace( natom_str );
      if (validInteger(natom_str)) {
        ptr = infile.NextLine(); // First atom
        if (ptr != 0) {
          char resnum[6], resname[6], atname[6], atnum[6];
          float XYZ[3];
          nread = sscanf(ptr, "%5c%5c%5c%5c%f %f %f", resnum, resname,
                         atname, atnum, XYZ, XYZ+1, XYZ+2);
        }
      }
    }
  }
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

// Traj_Gro::GetBox()
/** \return true if box coords detected in line, false otherwise.
  * \param ucell Array of length 9 set with unit cell vectors, X{xyz} Y{xyz} Z{xyz}
  * \param line Line to scan.
  */
bool Traj_Gro::GetBox(double* ucell, const char* line) const {
  float fXYZ[9];
  int nbox = sscanf( line, "%f %f %f %f %f %f %f %f %f", fXYZ, fXYZ+1, fXYZ+2,
                     fXYZ+3, fXYZ+4, fXYZ+5, fXYZ+6, fXYZ+7, fXYZ+8);
  if (nbox == 3) {
    // Box lengths. Convert nm->Ang. Assume orthogonal FIXME OK?
    ucell[0] = (double)fXYZ[0] * Constants::NM_TO_ANG;
    ucell[1] = 0;
    ucell[2] = 0;

    ucell[3] = 0;
    ucell[4] = (double)fXYZ[1] * Constants::NM_TO_ANG;
    ucell[5] = 0;

    ucell[6] = 0;
    ucell[7] = 0;
    ucell[8] = (double)fXYZ[2] * Constants::NM_TO_ANG;

    return true;
  } else if (nbox == 9) {
    // Box vectors. Convert nm->Ang.
    for (int i = 0; i != 9; i++)
      ucell[i] = fXYZ[i] * Constants::NM_TO_ANG;
    return true;
  } // Otherwise assume no box.
  return false;
}

// Traj_Gro::openTrajin()
int Traj_Gro::openTrajin() {
  currentSet_ = 0;
  return file_.OpenFileRead( fname_ );
}

// Traj_Gro::setupTrajin()
int Traj_Gro::setupTrajin(FileName const& fnameIn, Topology* trajParm)
{
  float fXYZ[9];
  fname_ = fnameIn; // TODO SetupRead for BufferedLine
  // Open file for reading 
  if (file_.OpenFileRead( fname_ )) return TRAJIN_ERR;
  // Read the title. May contain time value, 't= <time>'
  const char* ptr = file_.Line();
  if (ptr == 0) {
    mprinterr("Error: Reading title.\n");
    return TRAJIN_ERR;
  }
  std::string title( ptr );
  RemoveTrailingWhitespace(title);
  if (debug_ > 0)
    mprintf("\tTitle: %s\n", title.c_str());
  bool hasTime = true;
  // TODO Is it OK to assume there will never be a negative time value?
  double timeVal = GetTimeValue( ptr );
  if (timeVal < 0.0) hasTime = false;
  if (debug_ > 0)
    mprintf("\tTimeval= %g HasTime= %i\n", timeVal, (int)hasTime);
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
  ptr = file_.Line();
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
  // Attempt to read box line
  // v1(x) v2(y) v3(z) [v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)]
  ptr = file_.Line();
  Box groBox;
  if (ptr != 0) {
    double ucell[9];
    if (GetBox(ucell, ptr))
      groBox.SetupFromUcell( ucell );
  }
  // Set trajectory information. No temperature info.
  SetCoordInfo( CoordinateInfo(groBox, hasV, false, hasTime) );
  SetTitle( title );
  // Check for multiple frames. If nothing was read above, 1 frame, no box. If
  // box info was read but nothing more can be read, 1 frame. Otherwise there
  // are more frames.
  bool hasMultipleFrames = false;
  if (ptr != 0) {
    if (!groBox.HasBox()) // Assume another title read.
      hasMultipleFrames = true;
    else {
      ptr = file_.Line(); // Read line after box info
      if (ptr != 0)
        hasMultipleFrames = true;
    }
  }
  // Set up some info for performing blank reads.
  linesToRead_ = natom_;
  if (groBox.HasBox())
    linesToRead_ += 1;
  int nframes = 1;
  if (hasMultipleFrames) {
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
      for (int i = 0; i != linesToRead_; i++)
        ptr = file_.Line();
      if (ptr == 0) break;
      nframes++;
      ptr = file_.Line(); // Next title or EOF
    }
  }
  file_.CloseFile();
  return nframes;
}

/** Convert Gromacs velocity units (nm / ps) to Amber units (Ang / (1/20.455)ps). */
const double Traj_Gro::GMX_VEL_TO_AMBER = Constants::NM_TO_ANG / Constants::AMBERTIME_TO_PS;

// Traj_Gro::readFrame()
int Traj_Gro::readFrame(int fnum, Frame& frm) {
  if (fnum < currentSet_) {
    file_.CloseFile();
    file_.OpenFileRead( fname_ );
    currentSet_ = 0;
  }
  // Position file
  const char* ptr;
  for (int set = currentSet_; set != fnum; set++) {
    ptr = file_.Line(); // Title
    ptr = file_.Line(); // Natoms
    for (int i = 0; i != linesToRead_; i++)
      ptr = file_.Line(); // Atom (and possibly box)
    if (ptr == 0) return 1;
  }
  // Read current frame
  ptr = file_.Line(); // Title
  if (ptr == 0) return 1;
  if (CoordInfo().HasTime())
    frm.SetTime( GetTimeValue(ptr) );
  ptr = file_.Line(); // Natoms TODO check?
  double* Xptr = frm.xAddress();
  if (CoordInfo().HasVel()) {
    double* Vptr = frm.vAddress();
    for (int i = 0; i != natom_; i++, Xptr += 3, Vptr += 3) {
      ptr = file_.Line(); // Atom
      sscanf(ptr, "%*5c%*5c%*5c%*5c%lf %lf %lf %lf %lf %lf",
             Xptr, Xptr+1, Xptr+2, Vptr, Vptr+1, Vptr+2);
      for (int n = 0; n != 3; n++) {
        Xptr[n] *= Constants::NM_TO_ANG;
        Vptr[n] *= GMX_VEL_TO_AMBER;
      }
    }
  } else {
    for (int i = 0; i != natom_; i++, Xptr += 3) {
      ptr = file_.Line(); // Atom
      sscanf(ptr, "%*5c%*5c%*5c%*5c%lf %lf %lf", Xptr, Xptr+1, Xptr+2);
      for (int n = 0; n != 3; n++)
        Xptr[n] *= Constants::NM_TO_ANG;
    }
  }
  // Box read
  if (CoordInfo().HasBox()) {
    ptr = file_.Line();
    double ucell[9];
    if (GetBox(ucell, ptr))
      frm.ModifyBox().AssignFromUcell( ucell );
  }

  ++currentSet_;
  return 0;
}

// Traj_Gro::Info()
void Traj_Gro::Info() {
  mprintf("is a GRO file");
  if (CoordInfo().HasTime()) mprintf(", with time");
  if (CoordInfo().HasVel()) mprintf(", with velocities");
  if (CoordInfo().HasBox()) mprintf(", with box info"); 
}
