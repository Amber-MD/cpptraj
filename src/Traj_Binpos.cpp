#include "Traj_Binpos.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Traj_Binpos::Traj_Binpos() :
  bpnatom_(0),
  bpnatom3_(0),
  frameSize_(0),
  bpbuffer_(0)
{}

Traj_Binpos::~Traj_Binpos() {
  if (bpbuffer_!=0) delete[] bpbuffer_;
}

bool Traj_Binpos::ID_TrajFormat() {
  unsigned char buffer[4];
  buffer[0] = ' ';
  buffer[1] = ' ';
  buffer[2] = ' ';
  buffer[3] = ' ';
  if ( OpenFile() ) return false;
  IO->Read(buffer, 4);
  CloseFile();
  // Check for the magic header of the Scripps binary format.
  if (buffer[0] == 'f' &&
      buffer[1] == 'x' &&
      buffer[2] == 'y' &&
      buffer[3] == 'z')
    return true;
  return false;
}

int Traj_Binpos::openTraj() {
  unsigned char buffer[4]; 

  switch (access_) {
    case READ:
      if (OpenFile()) return 1;
      // Read past magic header
      IO->Read(buffer, 4);
      break;
    case APPEND:
      mprinterr("Error: Append not supported for binpos files.\n");
      return 1;
      break;
    case WRITE: 
      if (OpenFile()) return 1;
      // Always write header
      buffer[0] = 'f';
      buffer[1] = 'x';
      buffer[2] = 'y';
      buffer[3] = 'z';
      IO->Write(buffer, 4);
      break; 
  }
  return 0;
}

void Traj_Binpos::closeTraj() {
  CloseFile();
}

int Traj_Binpos::setupTrajin(Topology* trajParm) {
  // Open - reads past the 4 byte header
  if (openTraj()) return 1;

  // Binpos file is 4 byte header followed by binpos frame records.
  // Each frame record is an int (# atoms) followed by 3*natom
  // floats (the xyz coords).
  // First assume binpos file has consistent # of atoms in each record
  // and try to determine # of frames.
  IO->Read(&bpnatom_, sizeof(int));
  if (bpnatom_ != trajParm->Natom()) {
    mprinterr("Error: # of atoms in binpos file frame 1 (%i) is not equal to\n",bpnatom_);
    mprinterr("Error: the # of atoms in associated parm %s (%i)\n",
              trajParm->c_str(), trajParm->Natom());
    return -1;
  }
  bpnatom3_ = bpnatom_ * 3;
  frameSize_ = (size_t)bpnatom3_ * sizeof(float);
  off_t framesize = (off_t)frameSize_ + sizeof(int);
  off_t filesize = file_size_;
  if ( compressType_ != NO_COMPRESSION)
    filesize = uncompressed_size_;
  filesize -= 4; // Subtract 4 byte header
  int Frames = (int)(filesize / framesize);
  if ( (filesize % framesize) == 0 ) {
    seekable_ = true;
  } else {
    mprintf("Warning: %s: Could not accurately predict # frames. This usually \n",
            BaseName());
    mprintf("         indicates a corrupted trajectory. Will attempt to read %i frames.\n",
            Frames);
    seekable_=false;
  }

  mprintf("\t%i atoms, framesize=%lu, filesize=%lu, #Frames=%i\n", 
          bpnatom_, framesize, filesize, Frames);

  // Allocate space to read in floats
  if (bpbuffer_!=0)
    delete[] bpbuffer_;
  bpbuffer_ = new float[ bpnatom3_ ];

  closeTraj();
  return Frames;
}

int Traj_Binpos::readFrame(int set, double* X, double* V, double* box, double* T) {
  int natoms;
  off_t offset;

  if (seekable_) {
    offset = (off_t) set;
    offset *= (frameSize_ + sizeof(int));
    offset += 4;
    IO->Seek(offset);
  }
  // Read past natom
  if (IO->Read(&natoms, sizeof(int))<1)
    return 1;
  // Sanity check
  if ( natoms != bpnatom_ ) {
    mprinterr("Error: Reading of binpos files with varying # of atoms is not supported.\n");
    return 1;
  }
  // Read coords
  IO->Read(bpbuffer_, frameSize_);
  // Convert float to double
  for (int i = 0; i < bpnatom3_; ++i)
    X[i] = (double)bpbuffer_[i];
  return 0;
}

int Traj_Binpos::setupTrajout(Topology *trajParm, int NframesToWrite) {
  bpnatom_ = trajParm->Natom();
  bpnatom3_ = bpnatom_ * 3;
  frameSize_ = (size_t)bpnatom3_ * sizeof(float);
  // Allocate space to write out floats
  if (bpbuffer_!=0)
    delete[] bpbuffer_;
  bpbuffer_ = new float[ bpnatom3_ ];
  if (hasBox_) 
    mprintf("Warning: BINPOS format does not support writing of box coordinates.\n");

  return 0;
}

int Traj_Binpos::writeFrame(int set, double* X, double* V, double* box, double T) {
  IO->Write( &bpnatom_, sizeof(int) );
  // Convert double to float
  for (int i = 0; i < bpnatom3_; ++i)
    bpbuffer_[i] = (float)X[i];
  if (IO->Write( bpbuffer_, frameSize_ )) return 1;
  return 0;
}

void Traj_Binpos::info() {
  mprintf("is a BINPOS file");
}
 
