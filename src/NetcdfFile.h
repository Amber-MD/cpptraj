#ifndef INC_NETCDFFILE_H
#define INC_NETCDFFILE_H
#include <string>
#include "CoordinateInfo.h"
class Frame;
/// The base interface to NetCDF trajectory files.
class NetcdfFile {
  public:
    /// For determining NetCDF trajectory file type
    enum NCTYPE { NC_AMBERTRAJ = 0, NC_AMBERRESTART, NC_AMBERENSEMBLE, NC_UNKNOWN };
    /// For determining netcdf4/hdf5 vs netcdf3. Sync with NcFmtTypeStr_
    enum NC_FMT_TYPE { NC_NOTNC = 0, NC_V3, NC_V4 };
    /// \return Type of specified file.
    static NCTYPE GetNetcdfConventions(NC_FMT_TYPE&, const char*);
#   ifndef BINTRAJ
    NetcdfFile() { }
#   else
    /// CONSTRUCTOR
    NetcdfFile();
    /// \return Coordinate info corresponding to current setup. TODO have in variable?
    CoordinateInfo NC_coordInfo() const;
    /// Open NetCDF file for reading.
    int NC_openRead(std::string const&);
    /// Open previously created NetCDF file for writing.
    int NC_openWrite(std::string const&);
    /// Create NetCDF reservoir.
    int NC_createReservoir(bool, double, int, int&, int&);
    /// Create NetCDF trajectory file of given type.
    int NC_create(NC_FMT_TYPE, std::string const&, NCTYPE, int,
                  CoordinateInfo const&, std::string const&, int);
    /// Create NetCDF (v3) trajectory file of given type.
    int NC_create(std::string const&, NCTYPE, int, 
                  CoordinateInfo const&, std::string const&, int);
    /// Close NetCDF file, do not reset dimension/variable IDs.
    void NC_close();
    /// \return Title of NetCDF file.
    std::string const& GetNcTitle() const { return nctitle_; }
    /// Set up NetCDF file for reading.
    int NC_setupRead(std::string const&, NCTYPE, int, bool, bool, int);
    /// Read - Remd Values
    int ReadRemdValues(Frame&);
    /// Write - Remd Values
    int WriteRemdValues(Frame const&);
#   ifdef MPI
#   ifdef HAS_PNETCDF
    int parallelWriteRemdValues(int, Frame const&);
#   endif
#   endif
    /// Convert given float array to double.
    inline void FloatToDouble(double*,const float*) const;
    /// Convert given double array to float.
    inline void DoubleToFloat(float*,const double*) const; 
    /// DEBUG - Write start and count arrays to STDOUT
    void DebugIndices() const;
    /// DEBUG - Write all variable IDs to STDOUT
    void DebugVIDs() const;

    inline int Ncid()      const { return ncid_;                }
    inline int Ncatom()    const { return ncatom_;              }
    inline int Ncatom3()   const { return ncatom3_;             }
    inline int Ncframe()   const { return ncframe_;             }
    inline int CoordVID()  const { return coordVID_;            }
    /// \return Variable ID of integer-compressed coordinates
    int CompressedPosVID() const { return compressedPosVID_; }
  protected: // TODO: Make all private
    /// Enumerated type for variable IDs. MUST MATCH vidDesc_.
    enum VidType { V_COORDS = 0, V_VEL,  V_FRC, V_TEMP, V_BOXL,
                   V_BOXA,       V_TIME, V_IND, V_RIDX, V_CIDX,
                   V_REMDVALS,
                   NVID };
#   ifdef HAS_HDF5
    /// Set all variables compression level.
    int SetCompression(int);
    /// Set default frame chunk size
    int SetFrameChunkSize(int);
    /// Create an integer compressed version of the trajectory
    int NC_createIntCompressed(int);
    /// Write an integer-compressed frame to the trajectory
    int NC_writeIntCompressed(Frame const&);
    /// Read an integer-compressed frame from the trajectory
    int NC_readIntCompressed(int, Frame& frmIn);
#   endif
#   ifdef MPI
    void Broadcast(Parallel::Comm const&);
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
    std::string nctitle_;
  private:
    /// Strings corresponding to NC_FMT_TYPE
    static const char* NcFmtTypeStr_[];
    static const char* ConventionsStr_[];

    /// Descriptions of VidType. MUST MATCH VidType (except NVID).
    static const char* vidDesc_[];
    /// Enumerated type for dimension IDs. MUST MATCH didDesc_.
    enum DidType { D_FRAME = 0, D_ATOM, D_SPATIAL,
                   NDID }; // TODO everything else
    /// Descriptions of DidType. MUST MATCH DidType (except NDID).
    static const char* didDesc_[];

    /// \return NetCDF trajectory type based on conventions.
    static NCTYPE GetNetcdfConventions(int);
    /// Check NetCDF file conventions version.
    void CheckConventionsVersion();
    /// Write array[ncatom3] compressed
    int NC_writeIntCompressed(const double*, int);
    /// Read array[ncatom3] compressed
    int NC_readIntCompressed(double*, int, int);

    int NC_defineIcompressedVar(const char*, VidType, int&);

    bool Has_pH() const;
    bool HasRedOx() const;
    bool HasForces()       const { return (frcVID_ != -1);      }
    bool HasVelocities()   const { return (velocityVID_ != -1); }
    bool HasCoords()       const { return (coordVID_ != -1);    }
    bool HasTemperatures() const;
    bool HasTimes()        const { return (timeVID_ != -1);     }

    /// Read - Set up frame dimension ID and number of frames.
    int SetupFrameDim();
    /// Read - Set up ensemble dimension ID and number of members.
    int SetupEnsembleDim();
    /// Read - Set up coordinates, velocities, forces, # atoms
    int SetupCoordsVelo(bool, bool);
    /// Read - Set up time variable if present
    int SetupTime();
    /// Read - Set up box information if present.
    int SetupBox();
    /// Read - Set up temperature information if present.
    void SetupTemperature();
    /// Read - Set up replica index info if present.
    int SetupMultiD();

    /// Set compression level for variable ID (HDF5 only)
    int NC_setDeflate(VidType, int) const;
    /// Set frame chunk size for variable ID (HDF5 only)
    int NC_setFrameChunkSize(VidType, int) const;
#   ifdef HAS_HDF5
    /// Set desired compression level for variable ID.
    int SetCompression(VidType, int);
    /// Calculate integer compression factor of 10 from given power
    int calcCompressFactor(int);
    /// Increase variable chunk sizes
    int NC_setVarDimChunkSizes(VidType, int, int, std::vector<int> const&, int, std::vector<size_t>&) const;
#   endif

    int NC_defineTemperature(int*, int);
    inline void SetRemDimDID(int, int*) const;

#   ifdef HAS_HDF5
    std::vector<int> deflateLevels_; ///< Compression levels for each VID
    int compressedPosVID_;           ///< Coordinates integer VID
    int compressedVelVID_;           ///< Velocities integer VID
    double compressedFac_;           ///< Compression factor 
    std::vector<int> itmp_;          ///< Temp space for converting to int
    int fchunkSize_;                 ///< Frame chunk size
#   endif

    std::vector<double> RemdValues_; ///< Hold remd values
    ReplicaDimArray remDimType_;     ///< Type of each dimension (multi-D).
    ReplicaDimArray remValType_;     ///< Type of each value (single or multi-D).
    // TODO audit the dimension IDs, may not need to be class vars.
    Box nc_box_;          ///< Hold box information
    NCTYPE myType_;       ///< Current file type.
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
    int RemdValuesVID_;   ///< Replica values variable ID.
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
