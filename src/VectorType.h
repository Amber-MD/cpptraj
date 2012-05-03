#ifndef INC_VECTORTYPE_H
#define INC_VECTORTYPE_H
#include <string>
#include "AtomMask.h"
#include "ModesInfo.h"
#include "ArgList.h"
#include "DataSet.h"
#include "Topology.h"
#include "Frame.h"
#include "Lapack_Diag.h"
/// Hold coordinate vector information
//NOTE: Adapted from PTRAJ transformVectorInfo
class VectorType : public DataSet {
  public:
    enum vectorMode {
      VECTOR_NOOP=0,    VECTOR_PRINCIPAL_X, VECTOR_PRINCIPAL_Y, VECTOR_PRINCIPAL_Z,
      VECTOR_DIPOLE,    VECTOR_BOX,         VECTOR_MASK,        VECTOR_IRED,
      VECTOR_CORRPLANE, VECTOR_CORR,        VECTOR_CORRIRED
    };

    VectorType();
    ~VectorType();

    bool operator==(const VectorType&);
    int AssignMaster(VectorType*);
    int ReadModesFromFile();
    int Init(ArgList&);
    int Allocate(int);
    void Info();
    //void Print();
    int Setup(Topology*);

    int Size();
    int Xmax();
    void WriteBuffer(CharBuffer&, int);
    int Width();

    int Action_CORR(Frame*);
    int Action_DIPOLE(Frame*, Topology*);
    int Action_PRINCIPAL(Frame*);
    int Action_MASK( Frame* );
    int Action_IRED( Frame* );
    int Action_BOX( Frame* );

    vectorMode Mode() { return mode_; }
    bool NoModeInfo() { return modinfo_==0; }
  private:
    //std::string filename_;
    int totalFrames_;
    int frame_;
    vectorMode mode_;
    AtomMask mask_;
    AtomMask mask2_;
    //std::vector<Vec3> C_;
    //std::vector<Vec3> V_;
    double *cx_; 
    double *cy_;
    double *cz_;
    double *vx_; 
    double *vy_;
    double *vz_;
    Lapack_Diag Principal_;

    VectorType* master_;    ///< If 0 this vector has master ModesInfo 
    ModesInfo* modinfo_;    ///< Eigenmode info for CORRIRED
    std::string modesfile_;
    int ibeg_;
    int iend_;
    int order_;
    int npair_;
    double avgcrd_[3];
    double rave_;
    double r3iave_;
    double r6iave_;
    double *cftmp_;
    double *p2cftmp_;
    double *rcftmp_;

    void leastSquaresPlane(int, double *, double *, double *, double *);
    void sphericalHarmonics(int, int, double*, double, double[2]);
};
#endif
