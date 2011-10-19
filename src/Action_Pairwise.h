#ifndef INC_ACTIONS_PAIRWISE_H
#define INC_ACTIONS_PAIRWISE_H
/// Class: Pairwise 
/// Action to compare pairs of atoms. Will function in two ways:
/// 1) Calculate Lennard-Jones and Coulomb energy for each frame. Save
///    the total energy each frame, also save the total energy on each
///    atom. 
#include "Action.h"
class Pairwise: public Action {
    AtomMask Mask0;
    AtomMask RefMask;
    AmberParm *RefParm;
    Frame *RefFrame;
    bool *skipv;
    int *natexidx;
    bool hasExclusion;
    double kes;
    DataSet *ds_vdw;
    DataSet *ds_elec;
    double ELJ, Eelec;
    int Ncomparison;
    std::vector<double> ref_evdw;
    double cut_evdw, cut_evdw1;
    std::vector<double> atom_evdw;
    std::vector<double> ref_eelec;
    double cut_eelec, cut_eelec1;
    std::vector<double> atom_eelec;
    bool isReference;

    int AllocateExclusion(AmberParm *);
    int NumInteractions(AtomMask *, AmberParm *);
    int SetupExclusion(AmberParm *, int);
    double Energy_LJ(AmberParm *, int, int, double, double *);
    double Energy_Coulomb(AmberParm *, int, int, double, double *);
    void WriteCutFrame(AmberParm *, AtomMask *, double *, Frame *, char*);
    void Energy(AtomMask*,Frame*,AmberParm*);
  public:
    Pairwise();
    ~Pairwise();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
