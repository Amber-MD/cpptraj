#include "TrajectoryIO.h"
#ifdef MPI
int TrajectoryIO::SyncTrajIO(Parallel::Comm const& commIn) {
  if (coordInfo_.SyncCoordInfo(commIn)) return 1;
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
