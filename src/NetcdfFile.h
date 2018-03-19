#ifndef INC_NETCDFFILE_H
#define INC_NETCDFFILE_H
#include <string>
#include "Frame.h"
/// The base interface to NetCDF trajectory files.
class NetcdfFile {
  public:
    /// For determining NetCDF trajectory file type
    enum NCTYPE { NC_UNKNOWN = 0, NC_AMBERTRAJ, NC_AMBERRESTART, NC_AMBERENSEMBLE };
    /// \return Type of given file.
    NCTYPE GetNetcdfConventions(const char*);
#   ifndef BINTRAJ
    NetcdfFile() { }
#   else 
    NetcdfFile();
    /// \return NetCDF trajectory type based on conventions.
    NCTYPE GetNetcdfConventions();
    /// Check NetCDF file conventions version.
    void CheckConventionsVersion();
    /// \return Coordinate info corresponding to current setup. TODO have in variable?
    CoordinateInfo NC_coordInfo() const;
    /// Open NetCDF file for reading.
    int NC_openRead(std::string const&);
    /// Open previously created NetCDF file for writing.
    int NC_openWrite(std::string const&);
    /// Create NetCDF reservoir.
    int NC_createReservoir(bool, double, int, int&, int&);
    /// Create NetCDF trajectory file of given type.
    int NC_create(std::string const&, NCTYPE, int, 
                  CoordinateInfo const&, std::string const&);
    /// Close NetCDF file, do not reset dimension/variable IDs.
    void NC_close();
    /// \return Title of NetCDF file.
    std::string GetNcTitle() const;
    int NC_setupRead(NCTYPE, int, bool, bool);
    /// Read - Set up frame dimension ID and number of frames.
    int SetupFrameDim();
    /// Read - Set up ensemble dimension ID and number of members.
    int SetupEnsembleDim();
    /// Read - Set up coordinates, velocities, forces, # atoms
    int SetupCoordsVelo(bool, bool);
    /// Read - Set up time variable if present
    int SetupTime();
    /// Read - Set up box information if present.
    int SetupBox(NCTYPE);
    /// Read - Set up temperature information if present.
    void SetupTemperature();
    /// Read - Set up replica index info if present.
    int SetupMultiD();
    /// Read - Remd Values
    int ReadRemdValues(Frame&);
    /// Write - Remd Values
    int WriteRemdValues(Frame const&);
    /// Convert given float array to double.
    inline void FloatToDouble(double*,const float*) const;
    /// Convert given double array to float.
    inline void DoubleToFloat(float*,const double*) const; 
    /// DEBUG - Write start and count arrays to STDOUT
    void WriteIndices() const;
    /// DEBUG - Write all variable IDs to STDOUT
    void WriteVIDs() const;
    /// DEBUG - Write general debug info to STDOUT
    void NetcdfDebug() const;

    inline int Ncid()      const { return ncid_;                }
    inline int Ncatom()    const { return ncatom_;              }
    inline int Ncatom3()   const { return ncatom3_;             }
    inline int Ncframe()   const { return ncframe_;             }
    inline int CoordVID()  const { return coordVID_;            }
    bool HasForces()       const { return (frcVID_ != -1);      }
    bool HasVelocities()   const { return (velocityVID_ != -1); }
    bool HasCoords()       const { return (coordVID_ != -1);    }
    bool HasTemperatures() const;
    bool HasTimes()        const { return (timeVID_ != -1);     }
  protected: // TODO: Make all private
#   ifdef MPI
    void Sync(Parallel::Comm const&);
#   endif
    size_t start_[4];    ///< Array starting indices
    size_t count_[4];    ///< Array counts
    int ncid_;           ///< NetCDF file ID
    int ncframe_;        ///< Total number of frames in file
    int TempVID_;        ///< Temperature variable ID.
    int coordVID_;       ///< Coordinates variable ID.
    int velocityVID_;    ///< Velocity variable ID.
    int frcVID_;         ///< Force variable ID.
    int cellAngleVID_;   ///< Box angles variable ID.
    int cellLengthVID_;  ///< Box lengths variable ID.
    int timeVID_;        ///< Time variable ID.
    // MultiD REMD
    int remd_dimension_; ///< Number of replica dimensions.
    int indicesVID_;     ///< Variable ID for replica indices.
    int repidxVID_;      ///< Variable ID for overall replica index.
    int crdidxVID_;      ///< Variable ID for overall coordinate index.
    // NC ensemble
    int ensembleSize_;
  private:
    static const char* ConventionsStr_[];

    bool Has_pH() const;
    bool HasRedOx() const;

    int NC_defineTemperature(int*, int);
    inline void SetRemDimDID(NCTYPE, int, int*) const;

    std::vector<double> RemdValues_; ///< Hold remd values
    ReplicaDimArray remDimType_;     ///< Type of each dimension (multi-D).
    ReplicaDimArray remValType_;     ///< Type of each value (single or multi-D).

    Box nc_box_;          ///< Hold box information
    int ncdebug_;
    int ensembleDID_;     ///< Ensemble dimenison ID
    int frameDID_;        ///< Frames dimension ID
    int atomDID_;         ///< Atom dimension ID
    int ncatom_;          ///< Number of atoms
    int ncatom3_;         ///< Number of coordinates (# atoms * 3)
    int spatialDID_;      ///< Spatial dimension ID (3)
    int labelDID_;        ///< Box angle labels dimension ID (alpha, beta, gamma)
    int cell_spatialDID_; ///< Box lengths dimension ID
    int cell_angularDID_; ///< Box angles dimension ID
    int spatialVID_;      ///< Spatial (x, y, z) variable ID
    int cell_spatialVID_; ///< Box lengths variable ID
    int cell_angularVID_; ///< Box angles variable ID
    int RemdValuesVID_;  ///< Replica values variable ID.
#   endif /* BINTRAJ */
};
#ifdef BINTRAJ
// ----- Inline Functions ------------------------------------------------------
/** Convert float coords to double coords
  * NOTE: natom3 needs to match up with size of Coord!
  */
void NetcdfFile::FloatToDouble(double* X, const float* Coord) const {
  for (int i=0; i < ncatom3_; ++i)
    X[i] = (double)Coord[i];
}

/** Convert double coords to float coords */
void NetcdfFile::DoubleToFloat(float* Coord, const double* X) const {
  for (int i=0; i < ncatom3_; ++i)
    Coord[i] = (float)X[i];
}
#endif
#endif
