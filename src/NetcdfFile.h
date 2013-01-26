#ifndef INC_NETCDFFILE_H
#define INC_NETCDFFILE_H
#include <string>
/// The base interface to NetCDF trajectory files.
class NetcdfFile {
  public:
    /// For determining netcdf file type
    enum NCTYPE { NC_UNKNOWN = 0, NC_AMBERTRAJ, NC_AMBERRESTART };
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
    int NC_create(std::string const&,NCTYPE,int,bool,bool,bool,bool,std::string const&);
    void NC_close();

    int SetupFrame();
    int SetupCoordsVelo();
    int SetupTime();
    int SetupBox(double*,NCTYPE);
    int SetupTemperature();
    int SetupMultiD();

    void FloatToDouble(double*,const float*);
    void DoubleToFloat(float*,const double*); 

    inline int Ncid()    { return ncid_;    }
    inline int Ncatom()  { return ncatom_;  }
    inline int Ncatom3() { return ncatom3_; }
    inline int Ncframe() { return ncframe_; }
    bool HasVelocities() { return (velocityVID_ != -1); } 

    inline void SetNcatom( int natomIn ) { ncatom_ = natomIn; }
  protected: // TODO: Make all private
    size_t start_[3];
    size_t count_[3];

    int ncid_;
    int ncframe_;
    int TempVID_;
    int coordVID_;
    int velocityVID_;
    int cellAngleVID_;
    int cellLengthVID_;
    int timeVID_;
    // MultiD REMD
    int remd_dimension_;
    int indicesVID_;

    bool checkNCerr(int);
  private:
    int ncdebug_;
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

    double velocityScale_;

    std::string GetAttrText(int, const char *);
    int GetDimInfo(const char *, int *);
#   endif
};
#endif
