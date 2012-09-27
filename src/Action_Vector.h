#ifndef INC_ACTION_VECTOR_H
#define INC_ACTION_VECTOR_H
#include "Action.h"
#include "DataSet_Vector.h"
// DEBUG
//#incl ude "CpptrajFile.h"
//#incl ude "PDBfile.h"
/// Hold coordinate vector information
//NOTE: Adapted from PTRAJ transformVectorInfo
class Action_Vector : public Action {
  public:
    Action_Vector();

    void print();

    static void sphericalHarmonics(int, int, const double*, double, double[2]);
  private:
    DataSet_Vector* Vec_;
    std::string filename_;
    int totalFrames_;
    int frame_;
    AtomMask mask_;
    AtomMask mask2_;
    // DEBUG
    //PDBfile PDB;
    //CpptrajFile debugpdb;
    //CpptrajFile debuginert;

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

    Vec3 leastSquaresPlane(int, double *, double *, double *);
};
#endif
