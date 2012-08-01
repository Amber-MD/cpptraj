#ifndef INC_VECTORTYPE_H
#define INC_VECTORTYPE_H
#include <string>
#include "Action.h"
#include "ModesInfo.h"
// DEBUG
//#incl ude "CpptrajFile.h"
//#incl ude "PDBfile.h"
/// Hold coordinate vector information
//NOTE: Adapted from PTRAJ transformVectorInfo
class VectorType : public DataSet, public Action {
  public:
    enum vectorMode {
      VECTOR_NOOP=0,    VECTOR_PRINCIPAL_X, VECTOR_PRINCIPAL_Y, VECTOR_PRINCIPAL_Z,
      VECTOR_DIPOLE,    VECTOR_BOX,         VECTOR_MASK,        VECTOR_IRED,
      VECTOR_CORRPLANE, VECTOR_CORR,        VECTOR_CORRIRED
    };

    VectorType();
    ~VectorType();

    void print();

    // DataSet functions
    int Size();
    int Xmax();
    void WriteBuffer(CpptrajFile&, int);
    int Width();

    vectorMode Mode() { return mode_; }
    int Frame()       { return frame_; }
    int Order()       { return order_; }
    ModesInfo* ModeInfo() { return modinfo_; }
    // NOTE: Should calcs involving Cftmp etc be done internally?
    // Next 3 only used for 'analyze timecorr'
    double* Cftmp() { return cftmp_; }
    double* P2cftmp() { return p2cftmp_; }
    double* Rcftmp() { return rcftmp_; }

    // Currently only used for matrix IRED
    double Dot(const VectorType& rhs) {
      return (vx_[0]*rhs.vx_[0] + vy_[0]*rhs.vy_[0] + vz_[0]*rhs.vz_[0]);
    }
    // Currently only used for 'analyze timecorr'
    void PrintAvgcrd(CpptrajFile&);
  private:
    static const char ModeString[][12];

    std::string filename_;
    int totalFrames_;
    int frame_;
    vectorMode mode_;
    AtomMask mask_;
    AtomMask mask2_;
    double *cx_; 
    double *cy_;
    double *cz_;
    double *vx_; 
    double *vy_;
    double *vz_;

    VectorType* master_;    ///< If 0 this vector has master ModesInfo 
    ModesInfo* modinfo_;    ///< Eigenmode info for CORRIRED
    std::string modesfile_;
    int ibeg_;
    int iend_;
    int order_;
    int npair_;
    // Below only used in analyze timecorr
    double avgcrd_[3];
    double rave_;
    double r3iave_;
    double r6iave_;
    double *cftmp_;
    double *p2cftmp_;
    double *rcftmp_;
    // DEBUG
    //PDBfile PDB;
    //CpptrajFile debugpdb;
    //CpptrajFile debuginert;

    bool operator==(const VectorType&);
    int AssignMaster(VectorType*);
    int ReadModesFromFile();
    int Allocate(int);
    void Info();

    // Action functions
    int init();
    int setup();
    int action();

    int Action_CORR();
    int Action_DIPOLE();
    int Action_PRINCIPAL();
    int Action_MASK(  );
    int Action_IRED(  );
    int Action_BOX(  );

    void leastSquaresPlane(int, double *, double *, double *, double *);
    void sphericalHarmonics(int, int, double*, double, double[2]);
};
#endif
