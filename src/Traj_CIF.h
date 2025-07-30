#ifndef INC_TRAJ_CIF_H
#define INC_TRAJ_CIF_H
#include "TrajectoryIO.h"
#include "CIFfile.h"
// Class: Traj_CIF
/// TrajecttoryIO class for reading coordinates from CIF files.
class Traj_CIF : public TrajectoryIO {
  public:
    Traj_CIF() : Natoms_(0), Nmodels_(0), Cartn_x_col_(0),
                 Cartn_y_col_(0), Cartn_z_col_(0) {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_CIF(); }
  private:
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(FileName const&, Topology*);
    int openTrajin();
    int readFrame(int,Frame&);
    void Info();
    void closeTraj() {}
    int processWriteArgs(ArgList&, DataSetList const&) { return 0; }
    int readVelocity(int, Frame&)  { return 1; }
    int readForce(int, Frame&)     { return 1; }
    int processReadArgs(ArgList&)  { return 0; }
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&,int, bool) { return 1; }
    int writeFrame(int,Frame const&)                           { return 1; } 

    /// Determine name of block containing coordinates.
    int determineCoordsBlock();

    CIFfile file_;
    Box boxInfo_;
    int Natoms_;
    int Nmodels_;
    int Cartn_x_col_;
    int Cartn_y_col_;
    int Cartn_z_col_;
    std::string blockName_; ///< Name of block containing the coordinate information
    std::string entryName_; ///< Name of block containing entry info
    std::string xstr_;      ///< Name of column with X coords
    std::string ystr_;      ///< Name of column with Y coords
    std::string zstr_;      ///< Name of column with Z coords
    std::string nstr_;      ///< Name of the column containing atom number
};
#endif
