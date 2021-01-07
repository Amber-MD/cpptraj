#ifndef INC_ACTION_RADIAL_H
#define INC_ACTION_RADIAL_H
#include "Action.h"
#include "ImageOption.h"
/// Calculate the radial distribution (pair correlation) function.
class Action_Radial: public Action {
  public:
    Action_Radial();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Radial(); }
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
#   ifdef _OPENMP
    void CombineRdfThreads();
#   endif
    typedef std::vector<unsigned long> Iarray;

#   ifdef _OPENMP
    bool threadsCombined_;           ///< True if CombineRdfThreads() has been called.
#   endif
    ImageOption imageOpt_;           ///< Used to decide if imaging should be used.
    Iarray RDF_;                     ///< Hold bin counts.
    std::vector<Iarray> rdf_thread_; ///< Hold bin count on each thread.
    AtomMask Mask1_;                 ///< Atoms to calculate RDF for.
    AtomMask Mask2_;                 ///< Optional mask to calc RDF to atoms in Mask1.
    AtomMask OuterMask_;             ///< Mask with the most atoms.
    AtomMask InnerMask_;             ///< Mask with the fewest atoms.
    typedef std::vector<AtomMask> Marray;
    Marray Sites1_;
    Marray Sites2_;
    enum RmodeType { NORMAL=0, NO_INTRAMOL, CENTER1, CENTER2, BYSITE };
    RmodeType rmode_;                ///< Type of calculation to perform.
    enum SmodeType { OFF = 0, BYRES, BYMOL };
    SmodeType siteMode1_;
    SmodeType siteMode2_;
    Topology* currentParm_;   ///< Current topology, needed for NO_INTERMOL
    int intramol_distances_;  ///< # of intra-molecular distances for NO_INTERMOL.
    bool useVolume_;          ///< If true normalize based on input volume.
    double volume_;           ///< Hold sum of volume for averaging.
    double maximum2_;         ///< Largest distance squared that can be binned.
    double spacing_;          ///< Bin spacing.
    double one_over_spacing_; ///< 1/spacing, used to avoid man division ops.
    int numBins_;             ///< The number of bins.
    int numthreads_;          ///< Number of threads.
    unsigned long numFrames_; ///< Number of frames for which RDF is calcd.
    double density_;          ///< Particle density (molecules/Ang^3).
    DataSet* Dset_;
    DataSet* intrdf_;
    DataSet* rawrdf_;
    int debug_;

    int SetupSiteArrayByAtom(Marray&, AtomMask const&) const;
    int SetupSiteArrayByRes(Marray&, Topology const&, AtomMask const&) const;
    int SetupSiteArrayByMol(Marray&, Topology const&, AtomMask const&) const;
};
#endif  
