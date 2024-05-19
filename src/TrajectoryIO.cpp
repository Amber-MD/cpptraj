#include "TrajectoryIO.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
TrajectoryIO::TrajectoryIO() :
  debug_(0),
  nwarn_cell_not_xaligned_(0),
  nwarn_cell_not_symmetric_(0)
{}

/** Increment count of frames where it is an issue that the
  * unit cell is not X-aligned. Print a warning only the 
  * first time.
  */
void TrajectoryIO::incrementXalignWarnCount(int set, const char* desc) {
  if (nwarn_cell_not_xaligned_ == 0 || debug_ > 0)
    mprintf("Warning: Set %i; unit cell is not X-aligned. Box cannot be properly stored as %s.\n", set+1, desc);
  nwarn_cell_not_xaligned_++;
}

/** Increment count of frames where it is an issue that the
  * unit cell is not symmetric. Print a warning only the 
  * first time.
  */
void TrajectoryIO::incrementSymmetricWarnCount(int set, const char* desc) {
  if (nwarn_cell_not_symmetric_ == 0 || debug_ > 0)
    mprintf("Warning: Set %i; unit cell is not symmetric. Box cannot be properly stored as %s.\n", set+1, desc);
  nwarn_cell_not_symmetric_++;
}

/** Print any warnings. Ideally called just before or just after closeTraj() */
void TrajectoryIO::PrintWarnings(std::string const& fname) {
  if (nwarn_cell_not_xaligned_ > 0) {
    mprintf("Warning: %s: unit cell was not X-aligned for %u frames.\n", fname.c_str(), nwarn_cell_not_xaligned_);
    mprintf("Warning: This trajectory format expects unit cells to be X-aligned,\n"
            "Warning: so the unit cell orientation may not be correct.\n");
  }
  if (nwarn_cell_not_symmetric_ > 0) {
    mprintf("Warning: %s: unit cell axes were not symmetric for %u frames.\n", fname.c_str(), nwarn_cell_not_symmetric_);
    mprintf("Warning: This trajectory format expects unit cell info to be symmetric,\n"
            "Warning: so the unit cell orientation may not be correct.\n");
  }
  nwarn_cell_not_xaligned_ = 0;
  nwarn_cell_not_symmetric_ = 0;
}

#ifdef MPI
/** Broadcast trajectory IO info from master. */
int TrajectoryIO::BroadcastTrajIO(Parallel::Comm const& commIn) {
  if (coordInfo_.BroadcastCoordInfo(commIn)) return 1;
  int tSize = (int)title_.size();
  commIn.MasterBcast( &tSize, 1, MPI_INT );
  if (commIn.Master())
    commIn.MasterBcast( (char*)title_.c_str(), tSize, MPI_CHAR );
  else {
    char* titleIn = new char[ tSize + 1 ];
    commIn.MasterBcast( titleIn, tSize, MPI_CHAR );
    titleIn[tSize] = '\0';
    SetTitle( titleIn );
    delete[] titleIn;
  }
  return 0;
}
#endif
