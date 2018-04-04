#ifndef INC_DATASET_COORDS_TRAJIN_H
#define INC_DATASET_COORDS_TRAJIN_H
#include "DataSet_Coords.h"
#include "TrajectoryIO.h"
/// Used to read frames from a single trajectoryon disk.
/** Unlike DataSet_Coords_TRJ this does not have to keep
  * track of multiple trajectories and is safe to use
  * threaded.
  * TODO This is really just a stopgap until DataSet_Coords_TRJ
  * is made thread-safe.
  */
class DataSet_Coords_Trajin : public DataSet_Coords {
  public:
    DataSet_Coords_Trajin();
    ~DataSet_Coords_Trajin();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Coords_Trajin(); }
    int AddSingleTrajin(std::string const&, ArgList&, Topology*);
    // ---- DataSet functions -------------------
    size_t Size() const                          { return nframes_; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    void Info() const;
    void Add( size_t, const void* )              { return;            }
    int Allocate(SizeArray const&)               { return 0;          }
    // ----- DataSet_Coords functions ------------
    /// DISABLED: Add a frame.
    void AddFrame(Frame const& fIn) { }
    /// DISABLED: Set CRD at position with frame.
    void SetCRD(int idx, Frame const& fIn) { }
    /// Get a frame at position.
    void GetFrame(int idx, Frame& fIn) { Traj_->readFrame(idx, fIn); }
    /// Get a frame at position corresponding to mask. FIXME NOT THREAD SAFE
    void GetFrame(int idx, Frame& fIn, AtomMask const& mIn) {
      Traj_->readFrame(idx, readFrame_);
      fIn.SetCoordinates(fIn, mIn);
    }
    /// Set topology and coordinate information.
    int CoordsSetup(Topology const&, CoordinateInfo const&);
  private:
    TrajectoryIO* Traj_;
    Frame readFrame_;
    int nframes_;
};
#endif
