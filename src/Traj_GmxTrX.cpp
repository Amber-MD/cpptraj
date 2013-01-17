#include "Traj_GmxTrX.h"
#include "CpptrajStdio.h"
#include "ByteRoutines.h"

// CONSTRUCTOR
Traj_GmxTrX::Traj_GmxTrX() :
  isBigEndian_(false),
  format_(TRR),
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
  step_(0),
  nre_(0),
  precision_(0),
  dt_(0.0),
  lambda_(0.0) 
{}

void Traj_GmxTrX::GmxInfo() {
  mprintf("------------------------------\n");
  Info();
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
  mprintf("\tstep= %i\n", step_);
  mprintf("\tnre= %i\n", nre_);
  mprintf("\tprecision= %i\n", precision_);
  mprintf("\tdt= %f\n", dt_);
  mprintf("\tlambda= %f\n", lambda_);
}

//const unsigned char Traj_GmxTrX::Magic_TRR_[4] = {201, 7, 0, 0};
//const unsigned char Traj_GmxTrX::Magic_TRJ_[4] = {203, 7, 0, 0};
const int Traj_GmxTrX::Magic_ = 1993;

bool Traj_GmxTrX::IsTRX(CpptrajFile& infile) {
  int magic;
  if ( infile.Read( &magic, 4 ) != 4 ) return 1;
  if (magic != Magic_) {
    // See if this is big endian
    endian_swap( &magic, 1 );
    if (magic != Magic_) 
      return false;
    else
      isBigEndian_ = true;
  } else
    isBigEndian_ = false;
  // TODO: At this point file is trX, but not sure how best to differentiate 
  // between TRR and TRJ. For now do it based on extension. Default TRR.
  if      (infile.Extension() == ".trr") format_ = TRR;
  else if (infile.Extension() == ".trj") format_ = TRJ;
  else format_ = TRR; 
  return true;
}

bool Traj_GmxTrX::ID_TrajFormat(CpptrajFile& infile) {
  // File must already be set up for read
  if (infile.OpenFile()) return false;
  bool istrx = IsTRX(infile);
  infile.CloseFile();
  return istrx;
}

void Traj_GmxTrX::closeTraj() {
  file_.CloseFile();
}

int Traj_GmxTrX::read_int( int& ival ) {
  // ASSUMING 4 byte integers
  if ( file_.Read( &ival, 4 ) != 4 ) return 1;
  if (isBigEndian_) endian_swap( &ival, 1 );
  return 0;
}

int Traj_GmxTrX::read_real( float& fval ) {
  double dval;
  switch (precision_) {
    case sizeof(float):
      if (file_.Read( &fval, precision_ ) != precision_) return 1;
      if (isBigEndian_) endian_swap( &fval, 1 );
      break;
    case sizeof(double):
      if (file_.Read( &dval, precision_ ) != precision_) return 1;
      if (isBigEndian_) endian_swap8( &dval, 1 );
      fval = (float)dval;
      break;
    default:
      return 1;
  }
  return 0;
}

std::string Traj_GmxTrX::read_string( ) {
  int size = 0;
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

/** Open trX trajectory and read past header info. */
int Traj_GmxTrX::openTrajin() {
  unsigned char buffer[4];
  if (file_.OpenFile()) return 1;
  // Read past magic byte
  if (file_.Read(buffer, 4) != 4) return 1;
  // Read version for TRR
  int version = 0;
  if (format_ != TRJ)
    read_int( version );
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
  if ( read_int( step_ ) ) return 1;
  if ( read_int( nre_ ) ) return 1;
  // Determine precision
  if (x_size_ > 0)
    precision_ = x_size_ / (natoms_ * 3);
  else if (v_size_ > 0)
    precision_ = v_size_ / (natoms_ * 3);
  else if (f_size_ > 0)
    precision_ = f_size_ / (natoms_ * 3);
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
  if ( read_real( dt_ ) ) return 1;
  if ( read_real( lambda_ ) ) return 1;
  GmxInfo(); // DEBUG
  return 0;
}

int Traj_GmxTrX::setupTrajin(std::string const& fname, Topology* trajParm)
{
  if (file_.SetupRead( fname, debug_ )) return TRAJIN_ERR;
  // Call openTrajin, which will open and read past header
  if ( openTrajin() ) return TRAJIN_ERR;
  // Warn if # atoms in parm does not match
  if (trajParm->Natom() != natoms_) {
    mprinterr("Error: # atoms in TRX file (%i) does not match # atoms in parm %s (%i)\n",
              natoms_, trajParm->c_str(), trajParm->Natom());
    return TRAJIN_ERR;
  }
  return TRAJIN_ERR;
}

int Traj_GmxTrX::setupTrajout(std::string const& fname, Topology* trajParm,
                                 int NframesToWrite, bool append)
{
  return 1;
}

int Traj_GmxTrX::readFrame(int set,double *X, double *V,double *box, double *T) {
  return 1;
}

int Traj_GmxTrX::writeFrame(int set, double *X, double *V,double *box, double T) {
  return 1;
}

void Traj_GmxTrX::Info() {
  mprintf("is a GROMACS");
   if (format_ == TRR)
    mprintf(" TRR file,");
  else
    mprintf(" TRJ file,");
  if (isBigEndian_) 
    mprintf(" big-endian");
  else
    mprintf(" little-endian");
  mprintf("\n");
}
