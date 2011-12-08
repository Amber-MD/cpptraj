#ifndef INC_ACTIONS_PAIRWISE_H
#define INC_ACTIONS_PAIRWISE_H
#include "Action.h"
// Class: Pairwise 
/// Action to calculate nonbonded energy between pairs of atoms.
/** Functions in two ways:
  * - Calculate Lennard-Jones and Coulomb energy for each frame. 
  *   Also calculate the cumulative LJ and Coulomb energy on each atom.
  * - Calculate the Lennard-Jones and Coulomb energy between each
  *   pair of atoms in a reference structure. Calculate the difference
  *   in each pair from frame to reference (d = Ref - Frame). 
  */
class Pairwise: public Action {
    AtomMask Mask0;                 ///< Calculate energy for atoms in mask
    AtomMask RefMask;               ///< Reference mask
    AmberParm *RefParm;             ///< Reference parm
    Frame *RefFrame;                ///< Reference coordinates
    bool *skipv;                    ///< Set by SetupExclusion, T for excluded (skipped) atoms
    int *natexidx;                  ///< Indices into excluded atom list NATEX for each atom
    bool hasExclusion;              ///< If true use the exclusion list
    double kes;                     ///< Electrostatic constant, 1.0 when using Amber units
    DataSet *ds_vdw;                ///< Evdw dataset
    DataSet *ds_elec;               ///< Eelec dataset
    double ELJ, Eelec;              ///< Total Evdw and Eelec over all atoms
    int N_ref_interactions;         ///< # of pairwise interactions in reference
    std::vector<double> ref_evdw;   ///< Evdw for each pair of atoms in reference
    double cut_evdw, cut_evdw1;     ///< Evdw cutoff ( Evdw < cutevdw1 && Evdw > cutevdw )
    std::vector<double> atom_evdw;  ///< Cumulative Evdw on each atom
    std::vector<double> ref_eelec;  ///< Eelec for each pair of atoms in reference
    double cut_eelec, cut_eelec1;   ///< Eelec cutoff ( Eelec < cuteelec1 && Eelec > cuteelec )
    std::vector<double> atom_eelec; ///< Cumulative Eelec on each atom
    char *cutout;                   ///< Mol2 file prefix for atoms satisfying cutoffs
    CpptrajFile Eout;               ///< Output file for atom energies.

    int AllocateExclusion(AmberParm *);
    int NumInteractions(AtomMask *, AmberParm *);
    int SetupExclusion(AmberParm *, int);
    double Energy_LJ(AmberParm *, int, int, double, double *);
    double Energy_Coulomb(AmberParm *, int, int, double, double *);
    int WriteCutFrame(AmberParm *, AtomMask *, double *, Frame *, char*);
    void RefEnergy();
    int Energy(AtomMask*,Frame*,AmberParm*);
  public:
    Pairwise();
    ~Pairwise();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
