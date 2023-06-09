#ifndef INC_ANALYSIS_CONSTANTPHSTATS_H
#define INC_ANALYSIS_CONSTANTPHSTATS_H
#include <map>
#include "Analysis.h"
#include "DataSet_pH.h"
/// <Enter description of Analysis_ConstantPHStats here>
class Analysis_ConstantPHStats : public Analysis {
  public:
    Analysis_ConstantPHStats() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_ConstantPHStats(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    DataFile* fracPlotOut_;
    std::string dsname_;
    int debug_;
    bool createFracPlot_;
    bool useFracProtonated_;
    DataSetList inputSets_;
    DataSetList* masterDSL_;
#   ifdef MPI
    Parallel::Comm masterComm_; ///< Communicator containing only ensemble masters.
    /// \return True for parallel (but will only use ensemble masters) FIXME use all procs?
    bool IsParallel() const { return true; }
#   endif

    /// Hold all DataSets for a residue
    struct resStatData {
      DataSet* nTrans_;
      DataSet* nProt_;
      DataSet* totProt_;
      DataSet* fracProt_;
      DataSet* fracDeprot_;
    };
    /// Hold various protonation stats for input DataSet (single residue at single pH)
    class ResStat {
      public:
        /// CONSTRUCTOR
        ResStat() : ds_(0), n_transitions_(0), n_prot_(0), tot_prot_(0), nframes_(0) {}
        /// CONSTRUCTOR : pH DataSet
        ResStat(DataSet_pH* ds) :
          ds_(ds),
          n_transitions_(0), n_prot_(0), tot_prot_(0), nframes_(0) {}
        /// CONSTRUCTOR : pH DataSet, initial state
        ResStat(DataSet_pH* ds, unsigned int init_state) :
          ds_(ds),
          n_transitions_(0),
          n_prot_((unsigned int)ds->Res().IsProtonated(init_state)),
          tot_prot_(ds->Res().Nprotons(init_state)),
          nframes_(1)
        {}

        DataSet_pH* ds_;    ///< Pointer to the associated input DataSet.
        unsigned int n_transitions_; ///< Protonated -> Deprotonated or vice-versa
        unsigned int n_prot_;        ///< # states protonated
        unsigned int tot_prot_;      ///< Total proton count
        unsigned int nframes_;       ///< Total number of frames
    };
    typedef std::vector<ResStat> Rarray;

    /// Sort ResStat by pH, then residue number.
    struct ph_num_sort {
      inline bool operator()(ResStat const& r0, ResStat const& r1) const {
        if (r0.ds_->Solvent_pH() == r1.ds_->Solvent_pH())
          return (r0.ds_->Res().Num() < r1.ds_->Res().Num());
        else
          return (r0.ds_->Solvent_pH() < r1.ds_->Solvent_pH());
      }
    };

    /// Sort ResStat by residue number, then pH.
    struct num_ph_sort {
      inline bool operator()(ResStat const& r0, ResStat const& r1) const {
        if (r0.ds_->Res().Num() == r1.ds_->Res().Num())
          return (r0.ds_->Solvent_pH() < r1.ds_->Solvent_pH());
        else
          return (r0.ds_->Res().Num() < r1.ds_->Res().Num());
      }
    };

    
    Rarray Stats_;           ///< Hold residue protonation stats for each input DataSet
    CpptrajFile* statsOut_;
};
#endif
