#include <cmath>
#include "Traj_GmxTrX.h"
#include "CpptrajStdio.h"
#include "ByteRoutines.h"
#include "Constants.h" // PIOVER2

// CONSTRUCTOR
Traj_GmxTrX::Traj_GmxTrX() :
  swapBytes_(false),
  isBigEndian_(false),
  format_(TRR),
  dt_(1.0),
  ir_size_(0),
  e_size_(0),
  box_size_(0),
  vir_size_(0),
  pres_size_(0),
  top_size_(0),
  sym_size_(0),
  x_size_(0),
  v_size_(0),
  f_size_(0),
  natoms_(0),
  natom3_(0),
  step_(0),
  nre_(0),
  precision_(4),
  timestep_(0.0),
  lambda_(0.0),
  frameSize_(0),
  headerBytes_(0),
  arraySize_(0),
  farray_(0),
  darray_(0) 
{}

// DESTRUCTOR
Traj_GmxTrX::~Traj_GmxTrX() {
  if (farray_ != 0) delete[] farray_;
  if (darray_ != 0) delete[] darray_;
}

/** For debugging, print internal info. */
void Traj_GmxTrX::GmxInfo() {
  mprintf("------------------------------\nFile ");
  Info();
  mprintf("\n\tTitle= [%s]\n", Title().c_str());
  mprintf("\tir_size= %i\n", ir_size_);
  mprintf("\te_size= %i\n", e_size_);
  mprintf("\tbox_size= %i\n", box_size_);
  mprintf("\tvir_size= %i\n", vir_size_);
  mprintf("\tpres_size= %i\n", pres_size_);
  mprintf("\ttop_size= %i\n", top_size_);
  mprintf("\tsym_size= %i\n", sym_size_);
  mprintf("\tx_size= %i\n", x_size_);
  mprintf("\tv_size= %i\n", v_size_);
  mprintf("\tf_size= %i\n", f_size_);
  mprintf("\tnatoms= %i\n", natoms_);
  mprintf("\tnatom3= %i\n", natom3_);
  mprintf("\tstep= %i\n", step_);
  mprintf("\tnre= %i\n", nre_);
  mprintf("\tprecision= %i\n", precision_);
  mprintf("\tdt= %f\n", timestep_);
  mprintf("\tlambda= %f\n", lambda_);
}

//const unsigned char Traj_GmxTrX::Magic_TRR_[4] = {201, 7, 0, 0};
//const unsigned char Traj_GmxTrX::Magic_TRJ_[4] = {203, 7, 0, 0};
const int Traj_GmxTrX::Magic_ = 1993;

/** \return true if TRR/TRJ file. Determine endianness. */
bool Traj_GmxTrX::IsTRX(CpptrajFile& infile) {
  int magic;
  if ( infile.Read( &magic, 4 ) != 4 ) return 1;
  swapBytes_ = false;
  isBigEndian_ = false;
  if (magic != Magic_) {
    // See if this is big endian
    endian_swap( &magic, 1 );
    if (magic != Magic_) 
      return false;
    else {
      // Big-endian. If we are on little endian machine need to swap bytes.
      isBigEndian_ = true;
      if (!IsBigEndian()) swapBytes_ = true;
    }
  } else {
    // Little-endian (non-standard). If we are on big endian machine need to swap bytes.
    if (IsBigEndian()) swapBytes_ = true;
  }
  // TODO: At this point file is trX, but not sure how best to differentiate 
  // between TRR and TRJ. For now do it based on extension. Default TRR.
  if      (infile.Filename().Ext() == ".trr") format_ = TRR;
  else if (infile.Filename().Ext() == ".trj") format_ = TRJ;
  else format_ = TRR; 
  return true;
}

/** \return true if TRR/TRJ file. */
bool Traj_GmxTrX::ID_TrajFormat(CpptrajFile& infile) {
  // File must already be set up for read
  if (infile.OpenFile()) return false;
  bool istrx = IsTRX(infile);
  infile.CloseFile();
  return istrx;
}

// Traj_GmxTrX::closeTraj()
void Traj_GmxTrX::closeTraj() {
  file_.CloseFile();
}

/** Read 1 integer, swap bytes if necessary. */
int Traj_GmxTrX::read_int( int& ival ) {
  if ( file_.Read( &ival, 4 ) != 4 ) return 1;
  if (swapBytes_) endian_swap( &ival, 1 );
  return 0;
}

/** Write 1 integer, swap bytes if necessary. */
int Traj_GmxTrX::write_int( int ivalIn ) {
  int ival = ivalIn;
  // ASSUMING 4 byte integers
  if (swapBytes_) endian_swap( &ival, 1 );
  if ( file_.Write( &ival, 4 ) != 4) return 1;
  return 0;
}

/** Read 1 float/double based on precision, swap bytes if necessary. */
int Traj_GmxTrX::read_real( float& fval ) {
  double dval;
  switch (precision_) {
    case sizeof(float):
      if (file_.Read( &fval, precision_ ) != precision_) return 1;
      if (swapBytes_) endian_swap( &fval, 1 );
      break;
    case sizeof(double):
      if (file_.Read( &dval, precision_ ) != precision_) return 1;
      if (swapBytes_) endian_swap8( &dval, 1 );
      fval = (float)dval;
      break;
    default:
      return 1;
  }
  return 0;
}

/** Write 1 float/double based on precision, swap bytes if necessary. */
int Traj_GmxTrX::write_real( float fvalIn ) {
  float fval = fvalIn;
  double dval;
  switch (precision_) {
    case sizeof(float):
      if (swapBytes_) endian_swap( &fval, 1 );
      if (file_.Write( &fval, precision_ ) != precision_) return 1;
      break;
    case sizeof(double):
      dval = (double)fval;
      if (swapBytes_) endian_swap8( &dval, 1 );
      if (file_.Write( &dval, precision_ ) != precision_) return 1;
      break;
    default:
      return 1;
  }
  return 0;
}

/** Read an integer value which gives string size, then the following string
  * of that size.
  */
std::string Traj_GmxTrX::read_string( ) {
  int size = 0;
  const int BUF_SIZE = 128;
  char linebuffer_[BUF_SIZE];
  // Read string size
  if ( read_int( size ) ) return std::string();
  if ( size < BUF_SIZE ) {
    // Read entire string
    file_.Read( linebuffer_, size );
    linebuffer_[size] = '\0';
    return std::string(linebuffer_);
  } else {
    // String is larger than input buffer. Read into linebuffer until
    // entire string is read.
    std::string output;
    int chunksize = BUF_SIZE - 1;
    linebuffer_[chunksize] = '\0';
    int ntimes = size / chunksize;
    for (int i = 0; i < ntimes; i++) {
      file_.Read( linebuffer_, chunksize );
      output.append( linebuffer_ );
    }
    int leftover = size % chunksize;
    // Add any leftover
    if (leftover > 0) {
      file_.Read( linebuffer_, leftover );
      linebuffer_[leftover] = '\0';
      output.append( linebuffer_ );
    }
    return output;
  }
}

int Traj_GmxTrX::ReadTrxHeader() {
  int version = 0;
  // Read past magic byte
  if (file_.Read(&version, 4) != 4) return 1;
  // Read version for TRR
  if (format_ != TRJ)
    read_int( version );
  //mprintf("DEBUG: TRX Version= %i\n", version);
  // Read in title string
  SetTitle( read_string() );
  // Read in size data
  if ( read_int( ir_size_ ) ) return 1;
  if ( read_int( e_size_ ) ) return 1;
  if ( read_int( box_size_ ) ) return 1;
  if ( read_int( vir_size_ ) ) return 1;
  if ( read_int( pres_size_ ) ) return 1;
  if ( read_int( top_size_ ) ) return 1;
  if ( read_int( sym_size_ ) ) return 1;
  if ( read_int( x_size_ ) ) return 1;
  if ( read_int( v_size_ ) ) return 1;
  if ( read_int( f_size_ ) ) return 1;
  if ( read_int( natoms_ ) ) return 1;
  if (natoms_ < 1) {
    mprinterr("Error: No atoms detected in TRX trajectory.\n");
    return 1;
  }
  natom3_ = natoms_ * 3;
  if ( read_int( step_ ) ) return 1;
  if ( read_int( nre_ ) ) return 1;
  // Determine precision
  if (x_size_ > 0)
    precision_ = x_size_ / natom3_;
  else if (v_size_ > 0)
    precision_ = v_size_ / natom3_;
  else if (f_size_ > 0)
    precision_ = f_size_ / natom3_;
  else {
    mprinterr("Error: X/V/F sizes are 0 in TRX trajectory.\n");
    return 1;
  }
  if ( precision_ != sizeof(float) &&
       precision_ != sizeof(double) )
  {
    mprinterr("Error: TRX precision %i not recognized.\n", precision_);
    return 1;
  }
  // Read timestep and lambda
  if ( read_real( timestep_ ) ) return 1;
  if ( read_real( lambda_ ) ) return 1;
  return 0;
}

/** Open trX trajectory and read past header info. */
int Traj_GmxTrX::openTrajin() {
  if (file_.OpenFile()) return 1;
  return 0;
}

/** 
  * \param boxOut Double array of length 6 containing {X Y Z alpha beta gamma} 
  */
int Traj_GmxTrX::ReadBox(double* boxOut) {
  // xyz is an array of length 9 containing X{xyz} Y{xyz} Z{xyz}.
  double xyz[9];
  float f_boxIn[9];
  switch (precision_) {
    case sizeof(float):
      if (file_.Read( f_boxIn, box_size_ ) != box_size_) return 1;
      if (swapBytes_) endian_swap( f_boxIn, 9 );
      for (int i = 0; i < 9; ++i)
        xyz[i] = (double)f_boxIn[i];
      break;
    case sizeof(double):
      if (file_.Read( xyz, box_size_ ) != box_size_) return 1;
      if (swapBytes_) endian_swap8( xyz, 9 );
      break;
    default: return 1;
  }
  // Calculate box lengths
  // NOTE: GROMACS units are nm
  boxOut[0] = sqrt((xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2])) * 10.0;
  boxOut[1] = sqrt((xyz[3]*xyz[3] + xyz[4]*xyz[4] + xyz[5]*xyz[5])) * 10.0;
  boxOut[2] = sqrt((xyz[6]*xyz[6] + xyz[7]*xyz[7] + xyz[8]*xyz[8])) * 10.0;
  //mprintf("DEBUG:\tTRX Box Lengths: %f %f %f\n", boxOut[0], boxOut[1], boxOut[2]);
  if (boxOut[0] <= 0.0 || boxOut[1] <= 0.0 || boxOut[2] <= 0.0) {
    // Use zero-length box size and set angles to 90
    // TODO: This will cause box detection to fail in Trajin. Set to max(X,Y,Z)?
    boxOut[0] = boxOut[1] = boxOut[2] = 0.0;
    boxOut[3] = boxOut[4] = boxOut[5] = 90.0;
  } else {
    // Get angles between x+y(gamma), x+z(beta), and y+z(alpha)
    boxOut[5] = acos( (xyz[0]*xyz[3] + xyz[1]*xyz[4] + xyz[2]*xyz[5]) * 
                      100.0 / (boxOut[0]* boxOut[1]) ) * 90.0/Constants::PIOVER2;
    boxOut[4] = acos( (xyz[0]*xyz[6] + xyz[1]*xyz[7] + xyz[2]*xyz[8]) *
                      100.0 / (boxOut[0]* boxOut[2]) ) * 90.0/Constants::PIOVER2;
    boxOut[3] = acos( (xyz[3]*xyz[6] + xyz[4]*xyz[7] + xyz[5]*xyz[8]) *
                      100.0 / (boxOut[1]* boxOut[2]) ) * 90.0/Constants::PIOVER2;
  }
  //mprintf("DEBUG:\tTRX Box Angles: %f %f %f\n", boxOut[3], boxOut[4], boxOut[5]);
  return 0;
}

/** Prepare trajectory for reading. Determine number of frames. */
int Traj_GmxTrX::setupTrajin(FileName const& fname, Topology* trajParm)
{
  int nframes = 0;
  if (file_.SetupRead( fname, debug_ )) return TRAJIN_ERR;
  // Open and read in header
  if ( file_.OpenFile() ) return TRAJIN_ERR;
  ReadTrxHeader();
  if (debug_ > 0) GmxInfo(); // DEBUG
  // Warn if # atoms in parm does not match
  if (trajParm->Natom() != natoms_) {
    mprinterr("Error: # atoms in TRX file (%i) does not match # atoms in parm %s (%i)\n",
              natoms_, trajParm->c_str(), trajParm->Natom());
    return TRAJIN_ERR;
  }
  // If float precision, create temp array. Temp array not needed for double reads.
  if (precision_ == sizeof(float)) {
    if (farray_ != 0) delete[] farray_;
    farray_ = new float[ natom3_ ];
  } 
  // Attempt to determine # of frames in traj
  headerBytes_ = (size_t)file_.Tell();
  frameSize_ = headerBytes_ + (size_t)box_size_ + (size_t)vir_size_ + (size_t)pres_size_ +
                              (size_t)x_size_   + (size_t)v_size_ +   (size_t)f_size_;
                              //(size_t)ir_size_ + (size_t)e_size_ + (size_t)top_size_ + 
                              //(size_t)sym_size_;
  size_t file_size = (size_t)file_.UncompressedSize();
  if (file_size > 0) {
    nframes = (int)(file_size / frameSize_);
    if ( (file_size % frameSize_) != 0 ) {
      mprintf("Warning: %s: Number of frames in TRX file could not be accurately determined.\n"
              "Warning:   Will attempt to read %i frames.\n", file_.Filename().base(), nframes);
    }
  } else {
    mprintf("Warning: Uncompressed size could not be determined. This is normal for\n");
    mprintf("Warning: bzip2 files. Cannot check # of frames. Frames will be read until EOF.\n");
    nframes = TRAJIN_UNK;
  }
  // Load box info so that it can be checked.
  double box[6];
  box[0]=0.0; box[1]=0.0; box[2]=0.0; box[3]=0.0; box[4]=0.0; box[5]=0.0;
  if ( box_size_ > 0 ) {
    if ( ReadBox( box ) ) return TRAJIN_ERR;
  }
  // Set traj info - No time or temperature
  SetCoordInfo( CoordinateInfo(Box(box), (v_size_ > 0), false, false) );
  closeTraj();
  return nframes;
}

void Traj_GmxTrX::AllocateCoords() {
  // Allocate temp space for coords/velo
  if (farray_ != 0) {delete[] farray_; farray_ = 0;}
  if (darray_ != 0) {delete[] darray_; darray_ = 0;}
  arraySize_ = (size_t)natom3_;
  if (CoordInfo().HasVel()) arraySize_ *= 2;
  if (precision_ == sizeof(float)) 
    farray_ = new float[ arraySize_ ];
  else 
    darray_ = new double[ arraySize_ ];
}

void Traj_GmxTrX::WriteHelp() {
  mprintf("\tdt : Time step to multiply set #s by (default 1.0).\n");
}

int Traj_GmxTrX::processWriteArgs(ArgList& argIn) {
  dt_ = argIn.getKeyDouble( "dt", 1.0 );
  isBigEndian_ = true;
  if (!IsBigEndian()) swapBytes_ = true;
  // Prevent byte swapping (DEBUG)
  if (argIn.hasKey("noswap"))
    swapBytes_ = false;
  precision_ = sizeof(float);
  return 0;
}

int Traj_GmxTrX::setupTrajout(FileName const& fname, Topology* trajParm,
                              CoordinateInfo const& cInfoIn,
                              int NframesToWrite, bool append)
{
  if (!append) {
    SetCoordInfo( cInfoIn );
    natoms_ = trajParm->Natom();
    natom3_ = natoms_ * 3;
    // Default to TRR
    format_ = TRR;
    // Set up title
    if (Title().empty())
      SetTitle("Cpptraj generated TRR file.");
    // Set size defaults, box, velocity etc
    ir_size_ = 0;
    e_size_ = 0;
    if (CoordInfo().HasBox())
      box_size_ = precision_ * 9;
    else
      box_size_ = 0;
    vir_size_ = 0;
    pres_size_ = 0;
    top_size_ = 0;
    sym_size_ = 0;
    x_size_ = natom3_ * precision_;
    if (CoordInfo().HasVel())
      v_size_ = natom3_ * precision_;
    else
      v_size_ = 0;
    f_size_ = 0;
    step_ = 0;
    nre_ = 0;
    //dt_ = 0.0;
    lambda_ = 0.0;
    // Allocate temp space for coords/velo
    AllocateCoords();
    if (file_.SetupWrite( fname, debug_)) return 1;
    if (file_.OpenFile()) return 1;
  } else {
    int nframes = setupTrajin( fname, trajParm );
    if ( format_ == TRJ ) {
      mprinterr("Error: Only writes to TRR files supported.\n");
      return 1;
    }
    if ( nframes == TRAJIN_ERR ) return 1;
    mprintf("\tAppending to TRR file starting at frame %i\n", nframes);
    // Re-open for appending
    if (file_.SetupAppend( fname, debug_ )) return 1;
    if (file_.OpenFile()) return 1;
  }
  return 0;
}

/** Read array of size natom3 with set precision. Swap endianness if 
  * necessary. Since GROMACS units are nm, convert to Ang.
  */
int Traj_GmxTrX::ReadAtomVector( double* Dout, int size ) {
  switch (precision_) {
    case sizeof(float):
      if (file_.Read( farray_, size ) != size) return 1;
      if (swapBytes_) endian_swap(farray_, natom3_);
      for (int i = 0; i < natom3_; ++i)
        Dout[i] = (double)(farray_[i] * 10.0); // FIXME: Legit for velocities?
      break;
    case sizeof(double):
      if (file_.Read( Dout, size ) != size) return 1;
      if (swapBytes_) endian_swap8(Dout, natom3_);
      for (int i = 0; i < natom3_; ++i)
        Dout[i] *= 10.0; // FIXME: Legit for velocities?
      break;
    default: return 1;
  }
  return 0;
}

int Traj_GmxTrX::readFrame(int set, Frame& frameIn) {
  file_.Seek( (frameSize_ * set) + headerBytes_ );
  // Read box info
  if (box_size_ > 0) {
    if (ReadBox( frameIn.bAddress() )) return 1;
  }
  // Blank read past virial/pressure tensor
  file_.Seek( file_.Tell() + vir_size_ + pres_size_ );
  // Read coordinates
  if (x_size_ > 0) {
    if (ReadAtomVector(frameIn.xAddress(), x_size_)) {
      mprinterr("Error: Reading TRX coords frame %i\n", set+1);
      return 1;
    }
  }
  // Read velocities
  if (v_size_ > 0) {
    if (ReadAtomVector(frameIn.vAddress(), v_size_)) {
      mprinterr("Error: Reading TRX velocities frame %i\n", set+1);
      return 1;
    }
  }

  return 0;
}

int Traj_GmxTrX::readVelocity(int set, Frame& frameIn) {
  // Seek to frame and past box, virial, pressure, coords
  file_.Seek( (frameSize_ * set) + headerBytes_ + box_size_ + vir_size_ +
                                   pres_size_ + x_size_ );
  // Read velocities
  if (v_size_ > 0) {
    if (ReadAtomVector(frameIn.vAddress(), v_size_)) {
      mprinterr("Error: Reading TRX velocities frame %i\n", set+1);
      return 1;
    }
  }
  return 0;
}

int Traj_GmxTrX::writeFrame(int set, Frame const& frameOut) {
  int tsize;
  // Write header
  write_int( Magic_ );
  tsize = (int)Title().size() + 1;
  write_int( tsize );
  tsize = (int)Title().size();
  write_int( tsize );
  file_.Write( Title().c_str(), Title().size() );
  write_int( ir_size_ );
  write_int( e_size_ );
  write_int( box_size_ );
  write_int( vir_size_ );
  write_int( pres_size_ );
  write_int( top_size_ );
  write_int( sym_size_ );
  write_int( x_size_ );
  write_int( v_size_ );
  write_int( f_size_ );
  write_int( natoms_ );
  write_int( step_ );
  write_int( nre_ );
  float time = (float)(dt_ * (double)set);
  write_real( time ); // TODO: Write actual time
  write_real( lambda_ );
  // Write box
  // NOTE: GROMACS units are nm
  if (box_size_ > 0) {
    double ucell[9];
    double by = frameOut.BoxCrd().BoxY() * 0.1;
    double bz = frameOut.BoxCrd().BoxZ() * 0.1; 
    ucell[0] = frameOut.BoxCrd().BoxX() * 0.1;
    ucell[1] = 0.0;
    ucell[2] = 0.0;
    ucell[3] = by*cos(Constants::DEGRAD*frameOut.BoxCrd().Gamma());
    ucell[4] = by*sin(Constants::DEGRAD*frameOut.BoxCrd().Gamma());
    ucell[5] = 0.0;
    ucell[6] = bz*cos(Constants::DEGRAD*frameOut.BoxCrd().Beta());
    ucell[7] = (by*bz*cos(Constants::DEGRAD*frameOut.BoxCrd().Alpha()) - ucell[6]*ucell[3]) / 
                ucell[4];
    ucell[8] = sqrt(bz*bz - ucell[6]*ucell[6] - ucell[7]*ucell[7]);
    if (precision_ == sizeof(float)) {
      float f_ucell[9];
      for (int i = 0; i < 9; i++)
        f_ucell[i] = (float)ucell[i];
      if (swapBytes_) endian_swap( f_ucell, 9 );
      file_.Write( f_ucell, box_size_ );
    } else { // double
      if (swapBytes_) endian_swap8( ucell, 9 );
      file_.Write( ucell, box_size_ );
    }
  }
  // Write coords/velo
  // NOTE: GROMACS units are nm
  const double* Xptr = frameOut.xAddress();
  const double* Vptr = frameOut.vAddress();
  int ix = 0;
  if (precision_ == sizeof(float)) {
    for (; ix < natom3_; ix++)
      farray_[ix] = (float)(Xptr[ix] * 0.1);
    if (v_size_ > 0)
      for (int iv = 0; iv < natom3_; iv++, ix++)
        farray_[ix] = (float)(Vptr[iv] * 0.1);
    if (swapBytes_) endian_swap( farray_, arraySize_ );
    file_.Write( farray_, x_size_ + v_size_ );
  } else { // double
    for (; ix < natom3_; ix++)
      darray_[ix] = (Xptr[ix] * 0.1);
    if (v_size_ > 0)
      for (int iv = 0; iv < natom3_; iv++, ix++)
        darray_[ix] = (Vptr[iv] * 0.1);
    if (swapBytes_) endian_swap8( darray_, arraySize_ );
    file_.Write( darray_, x_size_ + v_size_ );
  }
  return 0;
}

void Traj_GmxTrX::Info() {
  mprintf("is a GROMACS");
   if (format_ == TRR)
    mprintf(" TRR file,");
  else
    mprintf(" TRJ file,");
  if (isBigEndian_) 
    mprintf(" big-endian,");
  else
    mprintf(" little-endian,");
  if (precision_ == sizeof(float))
    mprintf(" single precision");
  else if (precision_ == sizeof(double))
    mprintf(" double precision");
}
#ifdef MPI
// =============================================================================
int Traj_GmxTrX::parallelOpenTrajin(Parallel::Comm const& commIn) {
  mprinterr("Error: Parallel read not supported for GROMACS TRX.\n");
  return 1;
}

/** This assumes file has been previously set up with parallelSetupTrajout
  * and header has been written, so open append.
  */
int Traj_GmxTrX::parallelOpenTrajout(Parallel::Comm const& commIn) {
  return (file_.ParallelOpenFile( CpptrajFile::APPEND, commIn ));
}

/** First master performs all necessary setup, then sends info to all children.
  */
int Traj_GmxTrX::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                      CoordinateInfo const& cInfoIn,
                                      int NframesToWrite, bool append,
                                      Parallel::Comm const& commIn)
{
  int err = 0;
  // In parallel MUST know # of frames to write in order to correctly set size
  if (NframesToWrite < 1) {
    mprinterr("Error: # frames to write must be known for TRR output in parallel.\n");
    err = 1;
  } else if (commIn.Master()) {
    err = setupTrajout(fname, trajParm, cInfoIn, NframesToWrite, append);
    // Determine header size, (18 * 4) + titleSize TODO put in setupTrajout?
    headerBytes_ = (18 * 4) + Title().size();
    // Determine frame size
    frameSize_ = headerBytes_ + box_size_ + x_size_ + v_size_;
    // NOTE: setupTrajout leaves file open. Should this change?
    file_.CloseFile();
  }
  commIn.MasterBcast(&err, 1, MPI_INT);
  if (err != 0) return 1;
  // Synchronize info on non-master threads.
  SyncTrajIO( commIn );
  commIn.MasterBcast( &ir_size_, 1, MPI_INT );
  commIn.MasterBcast( &e_size_,  1, MPI_INT );
  commIn.MasterBcast( &box_size_, 1, MPI_INT );
  commIn.MasterBcast( &vir_size_, 1, MPI_INT );
  commIn.MasterBcast( &pres_size_, 1, MPI_INT );
  commIn.MasterBcast( &top_size_, 1, MPI_INT );
  commIn.MasterBcast( &sym_size_, 1, MPI_INT );
  commIn.MasterBcast( &x_size_, 1, MPI_INT );
  commIn.MasterBcast( &v_size_, 1, MPI_INT );
  commIn.MasterBcast( &f_size_, 1, MPI_INT );
  commIn.MasterBcast( &natoms_, 1, MPI_INT );
  commIn.MasterBcast( &natom3_, 1, MPI_INT );
  commIn.MasterBcast( &step_, 1, MPI_INT );
  commIn.MasterBcast( &nre_, 1, MPI_INT );
  commIn.MasterBcast( &precision_, 1, MPI_INT );
  commIn.MasterBcast( &timestep_, 1, MPI_DOUBLE );
  commIn.MasterBcast( &lambda_, 1, MPI_FLOAT );
  // NOTE: cast these to unsigned long long to avoid ambiguity since MPI doesnt have size_t
  unsigned long long buf[2];
  if (commIn.Master()) {
    buf[0] = (unsigned long long)frameSize_;
    buf[1] = (unsigned long long)headerBytes_;
    commIn.MasterBcast( buf, 2, MPI_UNSIGNED_LONG_LONG );
  } else {
    commIn.MasterBcast( buf, 2, MPI_UNSIGNED_LONG_LONG );
    frameSize_ = (size_t)buf[0];
    headerBytes_ = (size_t)buf[1];
    AllocateCoords(); // Should already be done on master
  }
  if (append)
    file_.SetupWrite( fname, debug_ );
  else
    file_.SetupAppend( fname, debug_ );
  if (debug_ > 0)
    rprintf("Gromacs TRR: parallel headerSize= %zu  frameSize= %zu\n", headerBytes_, frameSize_);

  return 0;
}

int Traj_GmxTrX::parallelReadFrame(int set, Frame& frameIn) { return 1; }

int Traj_GmxTrX::parallelWriteFrame(int set, Frame const& frameOut) {
  // Seek to given frame.
  file_.Seek( frameSize_ * set );
  return ( writeFrame(set, frameOut) );
}

void Traj_GmxTrX::parallelCloseTraj() { closeTraj(); }
#endif
