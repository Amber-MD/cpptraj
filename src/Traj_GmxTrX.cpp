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
  timestepPos_(0),
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
  if (isBigEndian_)
    mprintf("\tBig endian\n");
  else
    mprintf("\tLittle endian\n");
  if (swapBytes_)
    mprintf("\tSwapping bytes\n");
  else
    mprintf("\tNot swapping\n");
}

//const unsigned char Traj_GmxTrX::Magic_TRR_[4] = {201, 7, 0, 0};
//const unsigned char Traj_GmxTrX::Magic_TRJ_[4] = {203, 7, 0, 0};
const int Traj_GmxTrX::Magic_ = 1993;

const char* Traj_GmxTrX::Version_ = "GMX_trn_file";

/** Determine endianness and whether bytes need to be swapped. */
int Traj_GmxTrX::DetermineEndian(int magicIn) {
  int magic = magicIn;
  swapBytes_ = false;
  isBigEndian_ = false;
  if (magic != Magic_) {
    // See if this is big endian
    endian_swap( &magic, 1 );
    if (magic != Magic_) 
      return 1;
    else {
      // Big-endian. If we are on little endian machine need to swap bytes.
      isBigEndian_ = true;
      if (!IsBigEndian()) swapBytes_ = true;
    }
  } else {
    // Little-endian (non-standard). If we are on big endian machine need to swap bytes.
    if (IsBigEndian()) swapBytes_ = true;
  }
  return 0;
}

/** \return true if TRR/TRJ file. Determine endianness. */
bool Traj_GmxTrX::IsTRX(CpptrajFile& infile) {
  int magic;
  if ( infile.Read( &magic, 4 ) != 4 ) return false;
  if (DetermineEndian( magic )) return false;
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

// Traj_GmxTrX::ReadTrxHeader()
int Traj_GmxTrX::ReadTrxHeader(int& magic) {
  int version = 0;
  // Read past magic byte
  if (file_.Read(&magic, 4) != 4) return 1;
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
  // Save position just before timestep/lambda
  timestepPos_ = (size_t)file_.Tell();
  // Read timestep and lambda
  if ( read_real( timestep_ ) ) return 1;
  if ( read_real( lambda_ ) ) return 1;
  return 0;
}

/** Open trX trajectory. */
int Traj_GmxTrX::openTrajin() {
  if (file_.OpenFile()) return 1;
  return 0;
}

/** Read box information from current frame.
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
  boxOut[0] = sqrt((xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2])) * Constants::NM_TO_ANG;
  boxOut[1] = sqrt((xyz[3]*xyz[3] + xyz[4]*xyz[4] + xyz[5]*xyz[5])) * Constants::NM_TO_ANG;
  boxOut[2] = sqrt((xyz[6]*xyz[6] + xyz[7]*xyz[7] + xyz[8]*xyz[8])) * Constants::NM_TO_ANG;
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

// Traj_GmxTrX::AllocateCoords()
void Traj_GmxTrX::AllocateCoords() {
  // Allocate temp space for coords/velo/forces
  if (farray_ != 0) {delete[] farray_; farray_ = 0;}
  if (darray_ != 0) {delete[] darray_; darray_ = 0;}
  arraySize_ = (size_t)natom3_;
  if (v_size_ > 0) arraySize_ += (size_t)natom3_;
  if (f_size_ > 0) arraySize_ += (size_t)natom3_;
  if (debug_ > 0) {
    mprintf("DEBUG: Allocating array using precision %i\n", precision_);
    mprintf("DEBUG: arraySize is %zu\n", arraySize_);
  }
  if (precision_ == sizeof(float))
    farray_ = new float[ arraySize_ ];
  else 
    darray_ = new double[ arraySize_ ];
}

/** Prepare trajectory for reading. Determine number of frames. */
int Traj_GmxTrX::setupTrajin(FileName const& fname, Topology* trajParm)
{
  int nframes = 0;
  isBigEndian_ = true;
  if (!IsBigEndian()) swapBytes_ = true;
  if (file_.SetupRead( fname, debug_ )) return TRAJIN_ERR;
  // Open and read in header
  if ( file_.OpenFile() ) return TRAJIN_ERR;
  int magic;
  ReadTrxHeader(magic);
  if (DetermineEndian( magic )) {
    mprinterr("Error: File is not Gromacs TRR.\n");
    return TRAJIN_ERR;
  }
  if (debug_ > 0) GmxInfo(); // DEBUG
  // Warn if # atoms in parm does not match
  if (trajParm->Natom() != natoms_) {
    mprinterr("Error: # atoms in TRX file (%i) does not match # atoms in parm %s (%i)\n",
              natoms_, trajParm->c_str(), trajParm->Natom());
    return TRAJIN_ERR;
  }
  // Allocate temp coord arrays.
  AllocateCoords();
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
  // Box, coords, velocity, force, time 
  SetCoordInfo( CoordinateInfo(Box(box), true, (v_size_ > 0), (f_size_ > 0), true) );
  closeTraj();
  return nframes;
}

// Traj_GmxTrX::WriteHelp()
void Traj_GmxTrX::WriteHelp() {
  mprintf("\tdt : Time step to multiply set #s by (default 1.0). Ignored if time already present\n");
}

// Traj_GmxTrX::processWriteArgs()
int Traj_GmxTrX::processWriteArgs(ArgList& argIn) {
  dt_ = argIn.getKeyDouble( "dt", 1.0 );
  isBigEndian_ = true;
  if (!IsBigEndian()) swapBytes_ = true;
  // Prevent byte swapping (DEBUG)
  if (argIn.hasKey("noswap"))
    swapBytes_ = false;
  precision_ = sizeof(float);
  if (argIn.hasKey("double")) precision_ = sizeof(double);
  return 0;
}

// Traj_GmxTrX::setupTrajout()
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
      SetTitle(Version_);
    else
      mprintf("Warning: Using a custom title with TRR format may make the trajectory\n"
              "Warning:  incompatible with Gromacs analysis tools.\n");
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
    if (CoordInfo().HasForce())
      f_size_ = natom3_ * precision_;
    else
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

/** Convert Gromacs force units (kJ / mol * nm) to Amber units (kcal / mol * Ang) */
const double Traj_GmxTrX::GMX_FRC_TO_AMBER = Constants::ANG_TO_NM * Constants::J_TO_CAL;

/** Convert Amber force units to Gromacs */
const double Traj_GmxTrX::AMBER_FRC_TO_GMX = Constants::NM_TO_ANG * Constants::CAL_TO_J;

/** Convert Gromacs velocity units (nm / ps) to Amber units (Ang / (1/20.455)ps). */
const double Traj_GmxTrX::GMX_VEL_TO_AMBER = Constants::NM_TO_ANG / Constants::AMBERTIME_TO_PS;

/** Convert Amber velocity units to Gromacs */
const double Traj_GmxTrX::AMBER_VEL_TO_GMX = Constants::ANG_TO_NM * Constants::AMBERTIME_TO_PS;

// Traj_GmxTrX::readFrame()
int Traj_GmxTrX::readFrame(int set, Frame& frameIn) {
  file_.Seek( (frameSize_ * set) + timestepPos_ );
  // Read timestep and lambda
  if ( read_real( timestep_ ) ) return 1;
  if ( read_real( lambda_ ) ) return 1;
  frameIn.SetTime( timestep_ );
  // Read box info
  if (box_size_ > 0) {
    if (ReadBox( frameIn.bAddress() )) return 1;
  }
  // Blank read past virial/pressure tensor
  file_.Seek( file_.Tell() + vir_size_ + pres_size_ );
  // Read coords/velocities/forces
  int ix = 0;
  int total_size = x_size_ + v_size_ + f_size_;
  if (precision_ == sizeof(float)) {
    if (file_.Read( farray_, total_size ) != total_size) {
      mprinterr("Error: Could not read TRX frame %i\n", set+1);
      return 1;
    }
    if (swapBytes_) endian_swap(farray_, arraySize_);
    // Read coordinates
    if (x_size_ > 0) {
      double* Xptr = frameIn.xAddress();
      for (; ix != natom3_; ix++)
        Xptr[ix] = ((double)farray_[ix]) * Constants::NM_TO_ANG;
    }
    // Read velocities
    if (v_size_ > 0) {
      double* Vptr = frameIn.vAddress();
      for (int iv = 0; iv != natom3_; iv++, ix++)
        Vptr[iv] = ((double)farray_[ix]) * GMX_VEL_TO_AMBER;
    }
    // Read forces
    if (f_size_ > 0) {
      double* Fptr = frameIn.fAddress();
      for (int ir = 0; ir != natom3_; ir++, ix++)
        Fptr[ir] = ((double)farray_[ix]) * GMX_FRC_TO_AMBER;
    }
  } else if (precision_ == sizeof(double)) {
    if (file_.Read( darray_, total_size ) != total_size) {
      mprinterr("Error: Could not read TRX frame %i\n", set+1);
      return 1;
    }
    if (swapBytes_) endian_swap8(darray_, arraySize_);
    // Read coordinates
    if (x_size_ > 0) {
      double* Xptr = frameIn.xAddress();
      for (; ix != natom3_; ix++)
        Xptr[ix] = darray_[ix] * Constants::NM_TO_ANG;
    }
    // Read velocities
    if (v_size_ > 0) {
      double* Vptr = frameIn.vAddress();
      for (int iv = 0; iv != natom3_; iv++, ix++)
        Vptr[iv] = darray_[ix] * GMX_VEL_TO_AMBER;
    }
    // Read forces
    if (f_size_ > 0) {
      double* Fptr = frameIn.fAddress();
      for (int ir = 0; ir != natom3_; ir++, ix++)
        Fptr[ir] = darray_[ix] * GMX_FRC_TO_AMBER;
    }
  } else // SANITY CHECK
    mprinterr("Error: Unknown precision (%i)\n", precision_);

  return 0;
}

// Traj_GmxTrX::readVelocity()
int Traj_GmxTrX::readVelocity(int set, Frame& frameIn) {
  // Seek to frame and past box, virial, pressure, coords
  file_.Seek( (frameSize_ * set) + headerBytes_ + box_size_ + vir_size_ +
                                   pres_size_ + x_size_ );
  // Read velocities
  if (v_size_ > 0) {
    if (precision_ == sizeof(float)) {
      if (file_.Read( farray_, v_size_ ) != v_size_) {
        mprinterr("Error: Could not read velocities from TRX frame %i\n", set+1);
        return 1;
      }
      double* Vptr = frameIn.vAddress();
      for (int iv = 0; iv != natom3_; iv++)
        Vptr[iv] = ((double)farray_[iv]) * GMX_VEL_TO_AMBER;
    } else if (precision_ == sizeof(double)) {
      if (file_.Read( darray_, v_size_ ) != v_size_) {
        mprinterr("Error: Could not read velocities from TRX frame %i\n", set+1);
        return 1;
      }
      double* Vptr = frameIn.vAddress();
      for (int iv = 0; iv != natom3_; iv++)
        Vptr[iv] = darray_[iv] * GMX_VEL_TO_AMBER;
    }
  } else // SANITY
    mprintf("Warning: TRX file does not contain velocity information.\n");
  return 0;
}

// Traj_GmxTrX::readForce()
int Traj_GmxTrX::readForce(int set, Frame& frameIn) {
  // Seek to frame and past box, virial, pressure, coords, velo
  file_.Seek( (frameSize_ * set) + headerBytes_ + box_size_ + vir_size_ +
                                   pres_size_ + x_size_ + v_size_ );
  // Read forces 
  if (f_size_ > 0) {
    if (precision_ == sizeof(float)) {
      if (file_.Read( farray_, f_size_ ) != f_size_) {
        mprinterr("Error: Could not read forces from TRX frame %i\n", set+1);
        return 1;
      }
      double* Fptr = frameIn.fAddress();
      for (int ir = 0; ir != natom3_; ir++)
        Fptr[ir] = ((double)farray_[ir]) * GMX_FRC_TO_AMBER;
    } else if (precision_ == sizeof(double)) {
      if (file_.Read( darray_, f_size_ ) != f_size_) {
        mprinterr("Error: Could not read forces from TRX frame %i\n", set+1);
        return 1;
      }
      double* Fptr = frameIn.fAddress();
      for (int ir = 0; ir != natom3_; ir++)
        Fptr[ir] = darray_[ir] * GMX_FRC_TO_AMBER;
    }
  } else // SANITY
    mprintf("Warning: TRX file does not contain force information.\n");
  return 0;
}

// Traj_GmxTrX::writeFrame()
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
  float time;
  if (CoordInfo().HasTime()) 
    time = (float)frameOut.Time();
  else
    time = (float)(dt_ * (double)set);
  write_real( time );
  write_real( lambda_ );
  // Write box
  // NOTE: GROMACS units are nm
  if (box_size_ > 0) {
    Matrix_3x3 ucell = frameOut.BoxCrd().UnitCell( Constants::ANG_TO_NM );
    //mprintf("BoxX: %g %g %g BoxY: %g %g %g BoxZ: %g %g %g\n",
    //        ucell[0], ucell[1], ucell[2],
    //        ucell[3], ucell[4], ucell[5],
    //        ucell[6], ucell[7], ucell[8]);
    if (precision_ == sizeof(float)) {
      float f_ucell[9];
      for (int i = 0; i < 9; i++)
        f_ucell[i] = (float)ucell[i];
      if (swapBytes_) endian_swap( f_ucell, 9 );
      file_.Write( f_ucell, box_size_ );
    } else { // double
      if (swapBytes_) endian_swap8( ucell.Dptr(), 9 );
      file_.Write( ucell.Dptr(), box_size_ );
    }
  }
  // Write coords/velo/forces
  // NOTE: GROMACS units are nm
  const double* Xptr = frameOut.xAddress();
  const double* Vptr = frameOut.vAddress();
  const double* Fptr = frameOut.fAddress();
  int ix = 0;
  if (precision_ == sizeof(float)) {
    for (; ix < natom3_; ix++)
      farray_[ix] = (float)(Xptr[ix] * Constants::ANG_TO_NM);
    if (v_size_ > 0)
      for (int iv = 0; iv < natom3_; iv++, ix++)
        farray_[ix] = (float)(Vptr[iv] * AMBER_VEL_TO_GMX);
    if (f_size_ > 0)
      for (int ir = 0; ir < natom3_; ir++, ix++)
        farray_[ix] = (float)(Fptr[ir] * AMBER_FRC_TO_GMX);
    if (swapBytes_) endian_swap( farray_, arraySize_ );
    file_.Write( farray_, x_size_ + v_size_ + f_size_ );
  } else { // double
    for (; ix < natom3_; ix++)
      darray_[ix] = (Xptr[ix] * Constants::ANG_TO_NM);
    if (v_size_ > 0)
      for (int iv = 0; iv < natom3_; iv++, ix++)
        darray_[ix] = (Vptr[iv] * AMBER_VEL_TO_GMX);
    if (f_size_ > 0)
      for (int ir = 0; ir < natom3_; ir++, ix++)
        darray_[ix] = (Fptr[ir] * AMBER_FRC_TO_GMX);
    if (swapBytes_) endian_swap8( darray_, arraySize_ );
    file_.Write( darray_, x_size_ + v_size_ + f_size_ );
  }
  return 0;
}

// Traj_GmxTrX::Info()
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
  if (v_size_ > 0) mprintf(", velocities");
  if (f_size_ > 0) mprintf(", forces");
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
    frameSize_ = headerBytes_ + box_size_ + x_size_ + v_size_ + f_size_;
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
