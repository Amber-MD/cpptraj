#ifndef INC_NETCDFFILE_H
#define INC_NETCDFFILE_H
#include <string>
#include "CoordinateInfo.h"
/// The base interface to NetCDF trajectory files.
class NetcdfFile {
  public:
    /// For determining netcdf file type
    enum NCTYPE { NC_UNKNOWN = 0, NC_AMBERTRAJ, NC_AMBERRESTART, NC_AMBERENSEMBLE };
    NCTYPE GetNetcdfConventions(const char*);
#   ifndef BINTRAJ
    NetcdfFile() { }
#   else 
    NetcdfFile();

    void NetcdfDebug();
    std::string GetAttrText(const char *);
    NCTYPE GetNetcdfConventions();
    int NC_openRead(std::string const&);
    int NC_openWrite(std::string const&);
    int NC_createReservoir(bool, double, int, int&, int&);
    int NC_create(std::string const&, NCTYPE, int, 
                  CoordinateInfo const&, std::string const&);
    void NC_close();

    int SetupFrameDim();
    int SetupEnsembleDim();
    int SetupCoordsVelo(bool);
    int SetupTime();
    int SetupBox(double*,NCTYPE);
    int SetupTemperature();
    int SetupMultiD(ReplicaDimArray&);

    void FloatToDouble(double*,const float*);
    void DoubleToFloat(float*,const double*); 

    inline int Ncid()     const { return ncid_;    }
    inline int Ncatom()   const { return ncatom_;  }
    inline int Ncatom3()  const { return ncatom3_; }
    inline int Ncframe()  const { return ncframe_; }
    inline int CoordVID() const { return coordVID_; }
    bool HasForces()       { return (frcVID_ != -1); }
    bool HasVelocities()   { return (velocityVID_ != -1); }
    bool HasCoords()       { return (coordVID_ != -1);    }
    bool HasTemperatures() { return (TempVID_ != -1);     }
    bool HasTimes()        { return (timeVID_ != -1);     }
    inline void SetNcatom( int natomIn ) { ncatom_ = natomIn; }
    void WriteIndices() const;
    void WriteVIDs() const;
  protected: // TODO: Make all private
    size_t start_[4];
    size_t count_[4];

    int ncid_;
    int ncframe_;
    int TempVID_;             ///< Temperature variable ID.
    int coordVID_;            ///< Coordinates variable ID.
    int velocityVID_;         ///< Velocity variable ID.
    int frcVID_;              ///< Force variable ID.
    int cellAngleVID_;        ///< Box angles variable ID.
    int cellLengthVID_;       ///< Box lengths variable ID.
    int timeVID_;             ///< Time variable ID.
    // MultiD REMD
    int remd_dimension_;      ///< Number of replica dimensions.
    int indicesVID_;          ///< Variable ID for replica indices.
    int repidxVID_;           ///< Variable ID for overall replica index.
    int crdidxVID_;           ///< Variable ID for overall coordinate index.

    bool checkNCerr(int);
    void Sync();
  private:
    int ncdebug_;
    int ensembleDID_;
    int frameDID_;
    int atomDID_;
    int ncatom_;
    int ncatom3_;
    int spatialDID_;
    int labelDID_;
    int cell_spatialDID_;
    int cell_angularDID_;
    int spatialVID_;
    int cell_spatialVID_;
    int cell_angularVID_;

    std::string GetAttrText(int, const char *);
    int GetDimInfo(const char *, int *);
    int NC_defineTemperature(int*, int);
#   endif
};
#endif
