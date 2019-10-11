#ifndef NO_TNGFILE
#include <cmath>   // pow
#include <cstring> // strncmp
#include <cstdlib> // free
#include "Traj_TNG.h"
#include "CpptrajStdio.h"
#include "FileName.h"
#include "Topology.h"
#include "CpptrajFile.h"
#include "Constants.h"

/// CONSTRUCTOR
Traj_TNG::Traj_TNG() :
  values_(0),
  tngatoms_(0),
  tngframes_(-1),
  tngsets_(-1),
  current_frame_(-1),
  next_nblocks_(-1),
  next_blockIDs_(0),
  tngfac_(0),
  isOpen_(false)
{}

/// DESTRUCTOR
Traj_TNG::~Traj_TNG() {
  closeTraj();
  if (values_ != 0) free( values_ );
  if (next_blockIDs_ != 0) free( next_blockIDs_ );
}

/** Identify trajectory format. File should be setup for READ.
  * Note that as of version 1.8.2 the TNG library routines provide
  * no straightforward way to confirm that a file is indeed a TNG
  * file; providing a non-TNG file to the routines results in lots
  * of memory errors. Therefore, relying on the fact that the
  * 'GENERAL INFO' block should be present in all TNG files, and that
  * that string seems to always be in bytes 40-51.
  */
bool Traj_TNG::ID_TrajFormat(CpptrajFile& fileIn) {
  char tngheader[52];
  if (fileIn.OpenFile()) return false;
  if (fileIn.Read(&tngheader, 52) != 52) return false;
  fileIn.CloseFile();
  if (strncmp(tngheader+40, "GENERAL INFO", 12)==0) {
    mprintf("DEBUG: TNG file.\n");
    return true;
  }

  return false;
}

/** Print trajectory info to stdout. */
void Traj_TNG::Info() {
  mprintf("is a GROMACS TNG file");
}

/** Close file. */
void Traj_TNG::closeTraj() {
  mprintf("DEBUG: Calling closeTrajin() isOpen_=%1i\n", (int)isOpen_);
  if (isOpen_) {
    tng_util_trajectory_close(&traj_);
  }
  isOpen_ = false;
}

// -----------------------------------------------------------------------------
/** Open trajectory for reading. */
int Traj_TNG::openTrajin() {
  mprintf("DEBUG: Calling openTrajin() isOpen_=%1i\n", (int)isOpen_);
  if (isOpen_) 
    closeTraj();
  tng_function_status stat = tng_util_trajectory_open(filename_.full(), 'r', &traj_);
  if (stat != TNG_SUCCESS) {
    mprinterr("Error: Could not open TNG file '%s'\n", filename_.full());
    return TRAJIN_ERR;
  }
  mprintf("DEBUG: Successfully opened TNG file.\n");
  isOpen_ = true;
  current_frame_ = -1;

  return 0;
}

/** Read help */
void Traj_TNG::ReadHelp() {

}

/** Process read arguments. */
int Traj_TNG::processReadArgs(ArgList& argIn) {

  return 0;
}

/* Utility function for properly scaling coordinates according to the factor. */
void Traj_TNG::convertArray(double* out, float* in, unsigned int nvals) const {
  for (unsigned int i = 0; i != nvals; i++)
    out[i] = ((double)in[i]) * tngfac_;
}

/** \return 1 if no more blocks, -1 on error, 0 if ok.
  */
int Traj_TNG::getNextBlocks(int64_t &next_frame)
{
  tng_function_status stat = tng_util_trajectory_next_frame_present_data_blocks_find(
    traj_,
    current_frame_,
    blockIds_.size(),
    &blockIds_[0],
    &next_frame,
    &next_nblocks_,
    &next_blockIDs_);
  if (stat == TNG_CRITICAL) {
    return -1;
  } else if (stat == TNG_FAILURE) {
    return 1;
  }
  return 0;
}

static inline const char* DtypeStr(char typeIn) {
  switch (typeIn) {
    case TNG_INT_DATA : return "integer";
    case TNG_FLOAT_DATA : return "float";
    case TNG_DOUBLE_DATA : return "double";
    default : return "unknown";
  }
  return 0;
}

static inline const char* BtypeStr(int64_t typeIn) {
  switch (typeIn) {
    case TNG_TRAJ_BOX_SHAPE : return "box";
    case TNG_TRAJ_POSITIONS : return "positions";
    case TNG_TRAJ_VELOCITIES : return "velocities";
    case TNG_TRAJ_FORCES     : return "forces";
    case TNG_TRAJ_PARTIAL_CHARGES : return "partial charges";
    case TNG_TRAJ_FORMAL_CHARGES : return "formal charges";
    case TNG_TRAJ_B_FACTORS : return "B factors";
    case TNG_TRAJ_ANISOTROPIC_B_FACTORS : return "anisotropic B factors";
    case TNG_TRAJ_OCCUPANCY : return "occupancy";
    case TNG_TRAJ_GENERAL_COMMENTS : return "general comments";
    case TNG_TRAJ_MASSES : return "masses";
    default : return "unknown";
  }
  return 0;
}


int Traj_TNG::readValues(int64_t blockId, int64_t& next_frame, double& frameTime, char& datatype) {
  tng_function_status stat;
  int blockDependency;
  tng_data_block_dependency_get(traj_, blockId, &blockDependency);
  if (blockDependency & TNG_PARTICLE_DEPENDENT) {
    stat = tng_util_particle_data_next_frame_read( traj_,
                                                   blockId,
                                                   &values_,
                                                   &datatype,
                                                   &next_frame,
                                                   &frameTime );
  } else {
    stat = tng_util_non_particle_data_next_frame_read( traj_,
                                                       blockId,
                                                       &values_,
                                                       &datatype,
                                                       &next_frame,
                                                       &frameTime );
  }
  if (stat == TNG_CRITICAL) {
    mprinterr("Error: Could not read TNG block '%s'\n", BtypeStr(blockId));
    return -1;
  } else if (stat == TNG_FAILURE) {
    mprintf("Warning: Skipping TNG block '%s'\n", BtypeStr(blockId));
    return 1;
  }
  return 0;
}

/** Set up trajectory for reading.
  * \return Number of frames in trajectory.
  */
int Traj_TNG::setupTrajin(FileName const& fname, Topology* trajParm)
{
  filename_ = fname;
  // Open the trajectory
  if (openTrajin()) return TRAJIN_ERR;

  // Get number of atoms
  if (tng_num_particles_get(traj_, &tngatoms_) != TNG_SUCCESS) {
    mprinterr("Error: Could not get number of particles from TNG file.\n");
    return 1;
  }
  // Check number of atoms
  if (tngatoms_ != (int64_t)trajParm->Natom()) {
    mprinterr("Error: Number of atoms in TNG file (%li) does not match number\n"
              "Error:  of atoms in associated topology (%i)\n",
               tngatoms_, trajParm->Natom());
    return TRAJIN_ERR;
  }

  // Get number of frames
  tngframes_ = -1;
  if (tng_num_frames_get(traj_, &tngframes_) != TNG_SUCCESS) {
    mprinterr("Error: Could not get number of frames from TNG file.\n");
    return 1;
  }
  // Print number of frames
  mprintf("\tTNG file has %li frames.\n", tngframes_);

  // Get number of frame sets
  tngsets_ = -1;
  if (tng_num_frame_sets_get(traj_, &tngsets_) != TNG_SUCCESS) {
    mprinterr("Error: could not get number of frame sets from TNG file.\n");
    return 1;
  }
  mprintf("\tTNG file has %li frame sets.\n", tngsets_);
  int nframes = (int)tngsets_; // Could overflow here

  // Get the time per frame
  double tpf = 0.0;
  if (tng_time_per_frame_get( traj_, &tpf) != TNG_SUCCESS)
    tpf = 0.0;
  mprintf("\tTNG file time per frame: %g\n", tpf);

  // Get the exponential distance scaling factor
  int64_t tngexp;
  if (tng_distance_unit_exponential_get(traj_, &tngexp) != TNG_SUCCESS) {
    mprinterr("Error: Could not get distance scaling exponential from TNG.\n");
    return 1;
  }
  mprintf("\tTNG exponential: %li\n", tngexp);
  switch (tngexp) {
    case -9  :
      mprintf("\tTNG has units of nm\n");
      // Input is in nm. Convert to Angstroms.
      tngfac_ = 10.0;
      break;
    case -10 :
      mprintf("\tTNG has units of Angstrom\n");
      // Input is in Angstroms.
      tngfac_ = 1.0;
      break;
    default :
      // Convert to Angstroms.
      tngfac_ = pow(10.0, tngexp + 10); break;
  }
  // Print the scaling factor
  mprintf("\tTNG distance scaling factor: %g\n", tngfac_);

  // This will be used as temp space for reading in values from TNG
  if (values_ != 0) free( values_ ); 
  values_ = 0;
  if (next_blockIDs_ != 0) free( next_blockIDs_ );
  next_blockIDs_ = 0;

  // Set up blocks
  blockIds_.clear();
  blockIds_.push_back( TNG_TRAJ_BOX_SHAPE );
  blockIds_.push_back( TNG_TRAJ_POSITIONS );
  blockIds_.push_back( TNG_TRAJ_VELOCITIES );
  blockIds_.push_back( TNG_TRAJ_FORCES );
  blockIds_.push_back( TNG_GMX_LAMBDA );

  // Get block IDs in the first frame
  int64_t next_frame;
  if (getNextBlocks(next_frame) != 0) {
    mprinterr("Error: Unable to read blocks from first frame of TNG.\n");
    return TRAJIN_ERR;
  }
  bool hasVel = false;
  bool hasFrc = false;
  bool hasPos = false;
  Matrix_3x3 boxShape(0.0);
  for (int64_t idx = 0; idx < next_nblocks_; idx++)
  {
    int64_t blockId = next_blockIDs_[idx];
    mprintf("DEBUG: Block '%s' detected.\n", BtypeStr(blockId));
    if ( blockId == TNG_TRAJ_BOX_SHAPE ) {
      // Try to determine the box type by reading first frame values.
      double frameTime;
      char datatype;
      int err = readValues(blockId, next_frame, frameTime, datatype);
      if (err == -1) {
        mprinterr("Error: Critical error encountered reading block '%s'\n", BtypeStr(blockId));
        return TRAJIN_ERR;
      }
      // TODO: ok to only expect float data?
      if (datatype != TNG_FLOAT_DATA) {
        mprinterr("Error: TNG block '%s' data type is %s, expected float!\n", BtypeStr(blockId), DtypeStr(datatype));
        return TRAJIN_ERR;
      }
      convertArray(boxShape.Dptr(), (float*)values_, 9);
      boxShape.Print("boxShape");
    } else if ( blockId == TNG_TRAJ_POSITIONS ) {
      hasPos = true;
    } else if ( blockId == TNG_TRAJ_VELOCITIES ) {
      hasVel = true; 
    } else if ( blockId == TNG_TRAJ_FORCES ) {
      hasFrc = true; 
    }
  } // END loop over blocks

  // Set up coordinate info. Box, coord, vel, force, time
  SetCoordInfo( CoordinateInfo( Box(boxShape), hasPos, hasVel, hasFrc, (tpf > 0.0) ) );

  // Set up blocks that are actually there.
  blockIds_.clear();
  if (CoordInfo().HasBox())
    blockIds_.push_back( TNG_TRAJ_BOX_SHAPE );
  if (CoordInfo().HasCrd())
    blockIds_.push_back( TNG_TRAJ_POSITIONS );
  if (CoordInfo().HasVel())
    blockIds_.push_back( TNG_TRAJ_VELOCITIES );
  if (CoordInfo().HasForce())
    blockIds_.push_back( TNG_TRAJ_FORCES );

  closeTraj();
  return (int)nframes;
}

/** Read specified trajectory frame. */
int Traj_TNG::readFrame(int set, Frame& frameIn) {
  //int64_t numberOfAtoms = -1;
  //tng_num_particles_get(traj_, &numberOfAtoms); TODO could this change per frame?
  // Determine next frame with data
  int64_t next_frame;                   // Will get set to the next frame (MD) with data
  int64_t n_data_blocks_in_next_frame;  // Will get set to # blocks in next frame
  int64_t *block_ids_in_next_frame = 0; // Will hold blocks in next frame
  tng_function_status stat = tng_util_trajectory_next_frame_present_data_blocks_find(
    traj_,
    current_frame_,
    blockIds_.size(),
    &blockIds_[0],
    &next_frame,
    &n_data_blocks_in_next_frame,
    &block_ids_in_next_frame);
  if (stat == TNG_CRITICAL) {
    mprinterr("Error: could not get data blocks in next frame (set %i)\n", set+1);
    return 1;
  }
  if (stat == TNG_FAILURE) {
    mprintf("DEBUG: No more blocks.\n");
    return 1;
  }
  mprintf("DEBUG: Set %i next_frame %li nblocksnext %i\n", set+1, next_frame, n_data_blocks_in_next_frame);

  if (n_data_blocks_in_next_frame < 1) {
    mprinterr("Error: No data blocks in next frame (set %i, TNG frame %li)\n", set+1, next_frame);
    return 1;
  }

  // Process data blocks
  double frameTime;
  char datatype;
  for (int64_t idx = 0; idx < n_data_blocks_in_next_frame; idx++)
  {
    int64_t blockId = block_ids_in_next_frame[idx];
    int blockDependency;
    tng_data_block_dependency_get(traj_, blockId, &blockDependency);
    if (blockDependency & TNG_PARTICLE_DEPENDENT) {
      stat = tng_util_particle_data_next_frame_read( traj_,
                                                     blockId,
                                                     &values_,
                                                     &datatype,
                                                     &next_frame,
                                                     &frameTime );
    } else {
      stat = tng_util_non_particle_data_next_frame_read( traj_,
                                                         blockId,
                                                         &values_,
                                                         &datatype,
                                                         &next_frame,
                                                         &frameTime );
    }
    if (stat == TNG_CRITICAL) {
      mprinterr("Error: Could not read TNG block '%s'\n", BtypeStr(blockId));
    } else if (stat == TNG_FAILURE) {
      mprintf("Warning: Skipping TNG block '%s'\n", BtypeStr(blockId));
      continue;
    }
    mprintf("DEBUG: set %i frameTime %g block %s data %s\n", set+1, frameTime, BtypeStr(blockId), DtypeStr(datatype));
    // TODO: ok to only expect float data?
    if (datatype != TNG_FLOAT_DATA) {
      mprinterr("Error: TNG block '%s' data type is %s, expected float!\n", BtypeStr(blockId), DtypeStr(datatype));
      return 1;
    }
    // ----- Box -----------------------
    if ( blockId == TNG_TRAJ_BOX_SHAPE ) { // TODO switch?
      Matrix_3x3 boxShape(0.0);
      convertArray(boxShape.Dptr(), (float*)values_, 9);
      frameIn.SetBox( Box(boxShape) );
      mprintf("DEBUG: box set %i %g %g %g\n", set,
              frameIn.BoxCrd().BoxX(),
              frameIn.BoxCrd().BoxY(),
              frameIn.BoxCrd().BoxZ());
    // ----- Coords --------------------
    } else if ( blockId == TNG_TRAJ_POSITIONS ) {
      convertArray( frameIn.xAddress(), (float*)values_, tngatoms_*3 );
      const double* tmpXYZ = frameIn.XYZ(0);
      mprintf("DEBUG: positions set %i %g %g %g\n", set, tmpXYZ[0], tmpXYZ[1], tmpXYZ[2]);
    } else if ( blockId == TNG_TRAJ_VELOCITIES ) {
      convertArray( frameIn.vAddress(), (float*)values_, tngatoms_*3 );
    } else if ( blockId == TNG_TRAJ_FORCES ) {
      convertArray( frameIn.fAddress(), (float*)values_, tngatoms_*3 );
    }
  } // END loop over blocks in next frame
  if (block_ids_in_next_frame != 0) free(block_ids_in_next_frame);
  // TODO is it OK that frameTime is potentially set multiple times?
  frameIn.SetTime( frameTime / Constants::PICO );

  // Update current frame number
  current_frame_ = next_frame;

  return 0;
}

/** Read velocities from specified frame. */
int Traj_TNG::readVelocity(int set, Frame& frameIn) {

  return 0;
}

/** Read forces from specified frame. */
int Traj_TNG::readForce(int set, Frame& frameIn) {

  return 0;
}

// -----------------------------------------------------------------------------
/** Write help. */
void Traj_TNG::WriteHelp() {

}

/** Process write arguments. */
int Traj_TNG::processWriteArgs(ArgList& argIn, DataSetList const& DSLin) {

  return 0;
}

/** Set up trajectory for write. */
int Traj_TNG::setupTrajout(FileName const& fname, Topology* trajParm,
                                   CoordinateInfo const& cInfoIn, 
                                   int NframesToWrite, bool append)
{

  return 1;
}

/** Write specified trajectory frame. */
int Traj_TNG::writeFrame(int set, Frame const& frameOut) {

  return 0;
}

// =============================================================================
#ifdef MPI
/** Open trajectory for reading in parallel. */
int Traj_TNG::parallelOpenTrajin(Parallel::Comm const& commIn) {
  return 1;
}

/** Open trajectory for writing in parallel. */
int Traj_TNG::parallelOpenTrajout(Parallel::Comm const& commIn) {
  return 1;
}

/** Set up trajectory for write in parallel. */
int Traj_TNG::parallelSetupTrajout(FileName const& fname, Topology* trajParm,
                                           CoordinateInfo const& cInfoIn,
                                           int NframesToWrite, bool append,
                                           Parallel::Comm const& commIn)
{

  return 1;
}

/** Read frame in parallel. */
int Traj_TNG::parallelReadFrame(int set, Frame& frameIn) {

  return 1;
}

/** Write frame in parallel. */
int Traj_TNG::parallelWriteFrame(int set, Frame const& frameOut) {

  return 1;
}

/** Close trajectory in parallel. */
void Traj_TNG::parallelCloseTraj() {

}
#endif
#endif /* NO_TNGFILE */
