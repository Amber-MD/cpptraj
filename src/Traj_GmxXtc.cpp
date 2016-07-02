#include "Traj_GmxXtc.h"
#ifndef NO_XDRFILE
#include "CpptrajStdio.h"
#include "Constants.h"
#include <cstdio> // SEEK_SET, SEEK_CUR
#include "xdr_seek.h"

/// CONSTRUCTOR
Traj_GmxXtc::Traj_GmxXtc() : xd_(0), vec_(0), dt_(1.0), natoms_(0), prec_(1000) {}

/// DESTRUCTOR
Traj_GmxXtc::~Traj_GmxXtc() {
  closeTraj();
  if (vec_ != 0) delete[] vec_;
}

// Traj_GmxXtc::ID_TrajFormat()
bool Traj_GmxXtc::ID_TrajFormat(CpptrajFile& infile) {
  // See if we can read # atoms from xtc file.
  // NOTE: read_xtc_natoms *should* take const char* - boo
  if ( read_xtc_natoms( (char*)infile.Filename().full(), &natoms_ ) != exdrOK )
    return false;
  return (natoms_ > 0);
}

static inline void Next32ByteBoundary(int& offset) {
    offset += 3 - ((offset + 3) % 0x04);
}

// Traj_GmxXtc::setupTrajin()
int Traj_GmxXtc::setupTrajin(FileName const& fnameIn, Topology* trajParm)
{
  if (fnameIn.empty()) return TRAJIN_ERR;
  fname_ = fnameIn;
  // Read number of atoms
  if ( read_xtc_natoms( (char*)fname_.full(), &natoms_ ) != exdrOK ) {
    mprinterr("Error: Could not get number of atoms from XTC file.\n");
    return TRAJIN_ERR;
  }
  if (natoms_ != trajParm->Natom()) {
    mprinterr("Error: # atoms in XTC file (%i) does not match # atoms in parm %s (%i)\n",
              natoms_, trajParm->c_str(), trajParm->Natom());
    return TRAJIN_ERR;
  }
  // Allocate array for reading coords
  if (vec_ != 0) delete[] vec_;
  vec_ = new rvec[ natoms_ ];
  if (vec_ == 0) return TRAJIN_ERR;
  // Read one frame to determine box info
  if (openTrajin()) return TRAJIN_ERR;
  Frame tmp( natoms_ );
  // First frame offset is always zero
  frameOffsets_.push_back( 0 );
  if (readFrame(0, tmp)) return TRAJIN_ERR;
  // Determine total number of frames
  // FIXME in parallel every thread is doing this which is unnecessary.
  int nframes = TRAJIN_UNK;
  if ( natoms_ < 10 ) {
    // Small system, frame size is consistent
    CpptrajFile tmpfile;
    if (tmpfile.SetupRead(fname_, debug_)) return 1;
    off_t file_size = tmpfile.UncompressedSize();
    const off_t SMALL_HEADER_SIZE = 16 + (DIM*DIM*4) + 4; // 3*int+float, DIM^2 floats
    const off_t SMALL_BYTES_PER_ATOM = DIM * 4; // DIM floats
    off_t frame_size = SMALL_HEADER_SIZE + (SMALL_BYTES_PER_ATOM * natoms_);
    if ((file_size % frame_size) != 0) {
      mprinterr("Error: Could not determine number of frames in XTC file.\n");
      return TRAJIN_ERR;
    }
    nframes = (int)(file_size / frame_size);
    // Set up offsets
    frameOffsets_.reserve( nframes );
    for (off_t frm = 1; frm < nframes; frm++)
      frameOffsets_.push_back( frm * frame_size );
  } else {
    // Large system, use seek to determine # of frames
    const off_t LARGE_HEADER_SIZE = 44 + 8 + (DIM*DIM*4); // 11 int, 2 float, DIM^2 floats
    // Seek past header at beginning of file
    if (xdr_seek(xd_, LARGE_HEADER_SIZE, SEEK_SET) != 0) {
      mprinterr("Error: Could not seek to first frame in XTC.\n");
      return TRAJIN_ERR;
    }
    // Get first offset
    int offset;
    if (xdrfile_read_int( &offset, 1, xd_) == 0) {
      mprinterr("Error: Could not read first integer offset.\n");
      return TRAJIN_ERR;
    }
    // Advance offset to next frame boundary
    Next32ByteBoundary( offset );
    // Seek remaining frames
    nframes = 1;
    int err = 0;
    while (err == 0) {
      err = xdr_seek(xd_, (off_t)offset + LARGE_HEADER_SIZE, SEEK_CUR);
      if (err == 0) {
        // If no integer read it is likely the end of the file has been reached.
        if (xdrfile_read_int(&offset, 1, xd_) == 0) break;
        nframes++;
        // Store position
        frameOffsets_.push_back( xdr_tell(xd_) - 4 - LARGE_HEADER_SIZE );
        // Advance to next frame boundary
        Next32ByteBoundary( offset );
      }
    }
  }
  if (debug_ > 0)
    mprintf("DEBUG: %i frames, %zu offsets\n", nframes, frameOffsets_.size());
/*
  // Read number of frames
  int nframes = TRAJIN_UNK;
  unsigned long xtc_frames;
  if ( read_xtc_nframes( (char*)fname_.full(), &xtc_frames ) != exdrOK )
    mprintf("Warning: Could not determine # of frames in XTC file.\n");
  else
    nframes = (int)xtc_frames;
*/
  closeTraj();
  // No velocity, no temperature, yes time, no force
  SetCoordInfo( CoordinateInfo(ReplicaDimArray(), tmp.BoxCrd(), false, false, true, false) );
  return nframes;
}

// Traj_GmxXtc::setupTrajout()
int Traj_GmxXtc::setupTrajout(FileName const& fnameIn, Topology* trajParm,
                              CoordinateInfo const& cInfoIn,
                              int NframesToWrite, bool append)
{
  if (fnameIn.empty()) return 1;
  fname_ = fnameIn;
  if (!append) {
    SetCoordInfo( cInfoIn );
    natoms_ = trajParm->Natom();
    // Allocate array for writing coords
    if (vec_ != 0) delete[] vec_;
    vec_ = new rvec[ natoms_ ];
    if (vec_ == 0) return 1;
    // Open write
    xd_ = xdrfile_open(fname_.full(), "w");
    if (xd_ == 0) {
      mprinterr("Error: Could not open XTC file for write.\n");
      return 1;
    }
  } else {
    int nframes = setupTrajin( fname_, trajParm );
    if ( nframes == TRAJIN_ERR ) return 1;
    if ( nframes != TRAJIN_UNK )
      mprintf("\tAppending to XTC file starting at frame %i\n", nframes);
    // Re-open for append
    xd_ = xdrfile_open(fname_.full(), "a");
    if (xd_ == 0) {
      mprinterr("Error: Could not open XTC file for append.\n");
      return 1;
    }
  }
  return 0;
}

// Traj_GmxXtc::WriteHelp()
void Traj_GmxXtc::WriteHelp() {
  mprintf("\tdt : Time step to multiply set #s by (default 1.0). Ignored if time already present.\n");
}

// Traj_GmxXtc::processWriteArgs()
int Traj_GmxXtc::processWriteArgs(ArgList& argIn) {
  dt_ = argIn.getKeyDouble( "dt", 1.0 );
  return 0;
}

// Traj_GmxXtc::openTrajin()
int Traj_GmxXtc::openTrajin() {
  xd_ = xdrfile_open(fname_.full(), "r");
  if (xd_ == 0) {
    mprinterr("Error: Could not open XTC file for read.\n");
    return 1;
  }
  return 0;
}

// Traj_GmxXtc::closeTraj()
void Traj_GmxXtc::closeTraj() {
  if (xd_ != 0) xdrfile_close(xd_);
  xd_ = 0;
}

// Traj_GmxXtc::readFrame()
int Traj_GmxXtc::readFrame(int set, Frame& frameIn) {
  if (xdr_seek(xd_, frameOffsets_[set], SEEK_SET) != 0) {
    mprinterr("Error: Could not seek in XTC file, frame %i\n", set+1);
    return 1;
  }
  float time;
  int step;
  int result = read_xtc(xd_, natoms_, &step, &time, box_, vec_, &prec_);
  if (result != exdrOK) return 1;
  //mprintf("DEBUG: set %i step %i time %f\n", set, step, time);
  frameIn.SetTime( time );
  int idx = 0;
  for (int ix = 0; ix < natoms_; ix++)
    for (int kx = 0; kx < DIM; kx++)
      frameIn[idx++] = (double)vec_[ix][kx] * Constants::NM_TO_ANG;
  idx = 0;
  Matrix_3x3 ucell;
  for (int ii = 0; ii < DIM; ii++)
    for (int ij = 0; ij < DIM; ij++)
      ucell[idx++] = (double)box_[ii][ij] * Constants::NM_TO_ANG;
  frameIn.SetBox( Box(ucell) );
  return 0;
}

// Traj_GmxXtc::writeFrame()
int Traj_GmxXtc::writeFrame(int set, Frame const& frameOut) {
  float time;
  if (CoordInfo().HasTime())
    time = (float)frameOut.Time();
  else
    time = (float)(dt_ * (double)set);
  Matrix_3x3 Ucell = frameOut.BoxCrd().UnitCell( Constants::ANG_TO_NM );
  int idx = 0;
  for (int ii = 0; ii < DIM; ii++)
    for (int ij = 0; ij < DIM; ij++)
      box_[ii][ij] = (float)Ucell[idx++];
  idx = 0;
  for (int ix = 0; ix < natoms_; ix++)
    for (int kx = 0; kx < DIM; kx++)
      vec_[ix][kx] = (float)frameOut[idx++] * Constants::ANG_TO_NM;
  int result = write_xtc(xd_, natoms_, set, time, box_, vec_, prec_);
  if (result != 0) return 1;
  return 0;
}

// Traj_GmxXtc::Info()
void Traj_GmxXtc::Info() {
  mprintf("is a GROMACS XTC file");
}
#else /* NO_XDRFILE */
// =============================================================================
Traj_GmxXtc::Traj_GmxXtc() {}

Traj_GmxXtc::~Traj_GmxXtc() {}

bool Traj_GmxXtc::ID_TrajFormat(CpptrajFile&) { return false; }

int Traj_GmxXtc::setupTrajin(FileName const&, Topology*) { return 1; }

int Traj_GmxXtc::processWriteArgs(ArgList&) { return 0; }

int Traj_GmxXtc::setupTrajout(FileName const&,Topology*,CoordinateInfo const&,int,bool) {return 1;}

int Traj_GmxXtc::openTrajin() { return 1; }

void Traj_GmxXtc::closeTraj() { }

int Traj_GmxXtc::readFrame(int, Frame&) { return 1; }

int Traj_GmxXtc::writeFrame(int, Frame const&) { return 1; }

void Traj_GmxXtc::Info() { }

#endif /* NO_XDRFILE */
int Traj_GmxXtc::readVelocity(int set, Frame& frameIn) { return 1; }

int Traj_GmxXtc::readForce(int set, Frame& frameIn) { return 1; }

#ifdef MPI
// -----------------------------------------------------------------------------
int Traj_GmxXtc::parallelOpenTrajin(Parallel::Comm const& commIn) {
  mprinterr("Error: Parallel read not supported for GROMACS XTC.\n");
  return 1;
}

int Traj_GmxXtc::parallelOpenTrajout(Parallel::Comm const& commIn) { return 1; }

int Traj_GmxXtc::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                      CoordinateInfo const& cInfoIn,
                                      int NframesToWrite, bool append,
                                      Parallel::Comm const& commIn)
{
  return 1;
}

int Traj_GmxXtc::parallelReadFrame(int set, Frame& frameIn) { return 1; }

int Traj_GmxXtc::parallelWriteFrame(int set, Frame const& frameOut) { return 1; }

void Traj_GmxXtc::parallelCloseTraj() { }
#endif /* MPI */
