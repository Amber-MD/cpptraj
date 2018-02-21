#ifndef INC_REMDRESERVOIRNC_H
#define INC_REMDRESERVOIRNC_H
#include "NetcdfFile.h"
#include "Frame.h"
#include "FileName.h"
/// Use to write Amber NetCDF structure reservoirs with energy/temperature/bins
class RemdReservoirNC : private NetcdfFile {
  public:
    RemdReservoirNC() : eptotVID_(-1), binsVID_(-1) {}
    //void SetDebug(int d) { debug_ = d; } // TODO necessary?
    /// Initialize and open. Filename, title, coordinate info, # atoms, has bins, reservoir temp, seed
    int InitReservoir(FileName const&, std::string const&, CoordinateInfo const&, int, bool, double, int);
    /// Write structure, energy, and bin to reservoir
    int WriteReservoir(int, Frame const&, double, int);
    /// Close the reservoir
    void CloseReservoir() { NC_close(); }
  private:
    std::vector<float> Coord_; ///< For converting input coords to single precision
    int eptotVID_;             ///< Potential energy variable ID
    int binsVID_;              ///< Cluster bins variable ID
};
#endif
