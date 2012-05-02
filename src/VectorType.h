#ifndef INC_VECTORTYPE_H
#define INC_VECTORTYPE_H
#include <string>
#include "AtomMask.h"
#include "ModesInfo.h"
#include "ArgList.h"
#include "DataSet.h"
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
    int Init(ArgList&);
  private:
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

    bool master_;           ///< If true 
    ModesInfo modinfo_;
    int ibeg_;
    int iend_;
    int order_;
    int npair_;
    double *avgcrd_;
    double rave_;
    double r3iave_;
    double r6iave_;
    double *cftmp_;
    double *p2cftmp_;
    double *rcftmp_;
};
#endif
