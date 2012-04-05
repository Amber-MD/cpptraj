#ifndef INC_NETCDFFILE_H
#define INC_NETCDFFILE_H
#include <string>
class NetcdfFile {
  public:
    /// For determining netcdf file type
    enum NCTYPE { NC_UNKNOWN = 0, NC_AMBERTRAJ, NC_AMBERRESTART };
    NetcdfFile();
    NCTYPE GetNetcdfConventions();
#   ifdef BINTRAJ
    void NetcdfDebug();
    std::string GetAttrText(int, const char *);
    int GetDimInfo(const char *, int *);
    int NC_open(const char*);
    int NC_create(const char*);
    int SetupCoordinates();
    int SetupVelocity();
    int SetupTime();
    int SetupBox(double *);
    int SetupTemperature();
#   endif

  private:
#   ifdef BINTRAJ
    bool checkNCerr(int);
#   endif

    int ncdebug_;
    int ncid_;
    int atomDID_;
    int ncatom_;
    int ncatom3_;
    int coordVID_;
    int velocityVID_;
    int cellAngleVID_;
    int cellLengthVID_;
    int spatialDID_;
    int labelDID_;
    int cell_spatialDID_;
    int cell_angularDID_;
    int spatialVID_;
    int timeVID_;
    int cell_spatialVID_;
    int cell_angularVID_;
    int TempVID_;

    size_t start_[3];
    size_t count_[3];
};
#endif
