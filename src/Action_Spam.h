#ifndef INC_ACTION_SPAM_H
#define INC_ACTION_SPAM_H
#include "Action.h"
#include "ImagedAction.h"
#include "Vec3.h"
#include "Timer.h"
#include "PairList.h"
/**
SPAM is a water profiling technique developed by Guanglei Cui at
GlaxoSmithKline (GSK). The original implementation involved a set of specialized
Python scripts interacting with VMD (via the VolMap tool), numpy, NAMD (for the
SPAM energy calculations) and R (for the free energy calculation using a
specialized kernel density estimate). While that implementation demonstrated
proof of principle, simply re-ordering the trajectory for use with NAMD proved
to be a performance bottleneck because it was written in Python. SPAM was
rewritten from the ground up in cpptraj, significantly improving efficiency and
providing a simpler interface.

The original C++ implementation of SPAM in cpptraj was done by Jason Swails
while interning at GSK. This code was built as a patch on top of cpptraj v.12
and was rewritten by Jason Swails for the current cpptraj version.

 (C) 2012 - 2013
*/
class Action_Spam: public Action {
  public:
    Action_Spam();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Spam(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();
#   ifdef MPI
    int SyncAction();
    Parallel::Comm trajComm_;
#   endif

    typedef std::vector<int> Iarray;
    typedef std::vector<Vec3> Varray;
    typedef std::vector<Residue> Rarray;

    // ------------------- Functions -------------------
    int SetupParms(Topology const&);
    double Calculate_Energy(Frame const&, Residue const&);
    int Calc_G_Wat(DataSet*, unsigned int);
    // Custom Do- routines
    Action::RetType DoPureWater(int, Frame const&);
    Action::RetType DoSPAM(int, Frame&);

    typedef bool (Action_Spam::*FxnType)(Vec3, Vec3, double) const;
    bool inside_box(Vec3, Vec3, double) const;
    bool inside_sphere(Vec3, Vec3, double) const;
    inline double Ecalc(int, int, double) const;

    int debug_;
    FxnType Inside_;          ///< Function for determining if water is inside peak.
    ImagedAction image_;      ///< Imaging routines.
    PairList pairList_;       ///< Atom pair list (purewater_ only)
    std::vector<int> watidx_; ///< Hold water index for each atom (starting from 0).
    Matrix_3x3 ucell_;        ///< Unit cell matrix
    Matrix_3x3 recip_;        ///< Fractional matrix
    std::string solvname_;    ///< Name of the solvent residues
    double DG_BULK_;          ///< SPAM free energy of the bulk solvent
    double DH_BULK_;          ///< SPAM enthalpy of the bulk solvent
    double temperature_;      ///< Temperature at which SPAM simulation was run
    bool purewater_;          ///< True if running a pure water simulation to derive bulk properties
    bool reorder_;            ///< True if solvent should be reordered
    bool calcEnergy_;         ///< True if energy needs to be calculated.
    double cut2_;             ///< Non-bonded cutoff in Angstroms (squared)
    double onecut2_;          ///< 1 / cut2_
    double doublecut_;        ///< twice the cutoff (to test if boxes are big enough)
    CpptrajFile* infofile_;   ///< SPAM info file
    AtomMask mask_;           ///< Mask for selecting individual solvent residues
    Iarray resPeakNum_;       ///< Peak that each solvent residue is assigned to; -1 is unassigned
    std::string summaryfile_; ///< File containing the summary of all SPAM statistics
    double site_size_;        ///< Size of the water site. This is a full edge length or diameter
    std::vector<Iarray> peakFrameData_; ///< A list of all omitted frames for each peak
    Topology* CurrentParm_;             ///< Current topology (for NB params).
    std::vector<double> atom_charge_;   ///< Charges that have been converted to Amber units
    bool sphere_;                       ///< Is our site shape a sphere? If no, it's a box.
    DataSet* ds_dg_;                    ///< Hold final delta G values for each peak
    DataSet* ds_dh_;                    ///< Hold final delta H values for each peak
    DataSet* ds_ds_;                    ///< Hold final -T*S values for each peak
    std::vector<DataSet*> myDSL_;       ///< Hold energy data sets
    Varray peaks_;                      ///< List of each peak location
    Varray comlist_;                    ///< For given frame, each residue C.O.M. coords.
    Rarray solvent_residues_;           ///< List of each solvent residue
    int Nframes_;                       ///< Total number of frames
    bool overflow_;                     ///< True if cutoff overflowed our box coordinates
    // Timers
    Timer t_action_;
    Timer t_resCom_;
    Timer t_assign_;
    Timer t_occupy_;
    Timer t_energy_;
    Timer t_reordr_;
};

#endif
