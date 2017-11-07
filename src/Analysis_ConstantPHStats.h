#ifndef INC_ANALYSIS_CONSTANTPHSTATS_H
#define INC_ANALYSIS_CONSTANTPHSTATS_H
#include <map>
#include "Analysis.h"
#include "DataSet_PH.h"
/// <Enter description of Analysis_ConstantPHStats here>
class Analysis_ConstantPHStats : public Analysis {
  public:
    Analysis_ConstantPHStats() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_ConstantPHStats(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    int debug_;
    DataSetList inputSets_;

    /// Hold all DataSets for a residue
    struct resStatData {
      DataSet* nTrans_;
      DataSet* nProt_;
      DataSet* totProt_;
      DataSet* fracProt_;
      DataSet* fracDeprot_;
    };
    
    class ResStat {
      public:
        ResStat() : num_(-1), n_transitions_(0), n_prot_(0), tot_prot_(0) {}
        ResStat(DataSet_PH::Residue const& r, int init_state) :
          name_(r.Name()), num_(r.Num()), n_transitions_(0),
          n_prot_((int)r.IsProtonated(init_state)),
          tot_prot_(r.Nprotons(init_state))
        {}
/*        /// Sort by pH
        bool operator<(const ResStat& rhs) const {
          return (pH_ < rhs.pH_);
        }
        /// Equivalent if have same pH
        bool operator==(const ResStat&) const;*/

        NameType name_;
        int num_; // TODO remove?
//        float pH_;
        int n_transitions_;
        int n_prot_;
        int tot_prot_;
    };
    /// Map solvent pH to residue.
    typedef std::map<float,ResStat> PHresMap;
    /// Pair solvent pH to residue.
    typedef std::pair<float,ResStat> PHresPair;
    /// Map residue id to array of stats at different pH values.
    typedef std::map<int,PHresMap> StatMap;
    /// Pair residue id to array of stats
    typedef std::pair<int,PHresMap> StatPair;
};
#endif
