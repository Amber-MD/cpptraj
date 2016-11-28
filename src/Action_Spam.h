#ifndef INC_ACTION_SPAM_H
#define INC_ACTION_SPAM_H
#include "Action.h"
#include "ImagedAction.h"
#include "Vec3.h"
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
    // ------------------- Functions -------------------
    int SetupParms(Topology const&);
    double Calculate_Energy(Frame const&, Residue const&);
    void Calc_G_Wat(DataSet*, double&, double&, double&, double&, double&);
    // Custom Do- routines
    Action::RetType DoPureWater(int, Frame const&);
    Action::RetType DoSPAM(int, Frame&);

    typedef std::vector<int> Iarray;
    typedef std::vector<Vec3> Varray;

    typedef bool (Action_Spam::*FxnType)(Vec3, Vec3, double) const;
    bool inside_box(Vec3, Vec3, double) const;
    bool inside_sphere(Vec3, Vec3, double) const;

    FxnType Inside_;        ///< Function for determining if water is inside peak.
    ImagedAction image_;    ///< Imaging routines.
    std::string solvname_;  ///< Name of the solvent residues
    double bulk_;           ///< SPAM free energy of the bulk solvent
    bool purewater_;        ///< True if running a pure water simulation to derive bulk properties
    bool reorder_;          ///< True if solvent should be reordered
    bool calcEnergy_;       ///< True if energy needs to be calculated.
    double cut2_;           ///< Non-bonded cutoff in Angstroms (squared)
    double onecut2_;        ///< 1 / cut2_
    double doublecut_;      ///< twice the cutoff (to test if boxes are big enough)
    CpptrajFile* infofile_; ///< SPAM info file
    AtomMask mask_;         ///< Mask for selecting individual solvent residues
    std::string summaryfile_; ///< File containing the summary of all SPAM statistics
    double site_size_;        ///< Size of the water site. This is a full edge length or diameter
    /// A list of all omitted frames for each peak
    std::vector<Iarray> peakFrameData_;
    /// For a given frame, peak that each solvent residue is assigned to; -1 is unassigned
    Iarray resPeakNum_;
    /// Topology, for extracting necessary parameters for energy evaluations
    Topology* CurrentParm_;
    std::vector<double> atom_charge_; ///< List of charges that have been converted to Amber units
    bool sphere_;                     ///< Is our site shape a sphere? If no, it's a box.
    std::vector<DataSet*> myDSL_;     ///< Hold double precision data sets
    Varray peaks_;                    ///< List of each peak location
    Varray comlist_;                  ///< For given frame, each residue C.O.M. coords.
    /// List of the first atom and last atoms of each solvent residue
    std::vector<Residue> solvent_residues_;
    int Nframes_;                     ///< Total number of frames
    bool overflow_;                   ///< True if cutoff overflowed our box coordinates
    // DEBUG
    int set_counter_;
};

#endif
