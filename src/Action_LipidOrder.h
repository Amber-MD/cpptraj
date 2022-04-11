#ifndef INC_ACTION_LIPIDORDER_H
#define INC_ACTION_LIPIDORDER_H
#include "Action.h"
#include "CharMask.h"
/// Calculate lipid order parameters. 
class Action_LipidOrder : public Action {
  public:
    Action_LipidOrder();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_LipidOrder(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();

    static const unsigned int MAX_H_;

    enum AxisType { DX = 0, DY, DZ };

    /// Pair residue name and atom name; used to ID chains
    typedef std::pair<NameType, NameType> Npair;
    typedef std::vector<Npair> Narray;

    typedef std::vector<int> Iarray;
    typedef std::vector<double> Darray;

    /// Hold data for unique carbon type
    class CarbonData;

    /// Lipid chain data
    typedef std::vector<CarbonData> ChainType;

    /// Array of different chains
    typedef std::vector<ChainType> ChainArray;

    /// Hold information for a single lipid carbon site
    class CarbonSite;
    typedef std::vector<CarbonSite> Carray;

    /// Find existing chain type by res name/atom name or add new chain type.
    int FindChain(Npair const&);

    CharMask mask_;
    Narray Types_;      ///< Array of chain types (res name, atom name)
    ChainArray Chains_; ///< Hold unique chains
    Iarray Nchains_;    ///< How many of each chain type?
    Carray Sites_;      ///< Hold all carbon sites
    AxisType axis_;
    std::string dsname_;
    DataSetList* masterDSL_;
    DataFile* outfile_;
    int debug_;
    bool report_p2_;       ///< If true report raw <P2> instead of |<P2>|
#   ifdef MPI
    int SyncAction();
    Parallel::Comm trajComm_;
#   endif
    int nthreads_;
};

/// Hold data for carbon position in a chain.
class Action_LipidOrder::CarbonData {
  public:
    CarbonData();

    unsigned int Nvals()   const { return nvals_; }
    bool Init()            const { return init_;  }
    NameType const& Name() const { return name_;  }
    const char* name()     const { return *name_; }
    unsigned int NumH()    const { return nH_;    }
    double Avg(int, double&) const;

    void UpdateAngle(int i, double val) {
      sum_[i] += val;
      sum2_[i] += val * val;
    }
    void SetName( NameType const& n) {
      name_ = n;
      init_ = true;
    }
    void UpdateNvals()                    { nvals_++;                   }
    void SetNumH(unsigned int n)          { nH_ = n;                    }
#   ifdef _OPENMP
    /// Each thread needs its own scratch space
    void AllocateTempSpace(int n)         { tempScd_.assign( nH_*n, 0 ); }
    /// Sum scratch space back to main thread
    void ConsolidateTemp(int n) {
      for (int thread = 1; thread < n; thread++) {
        int offset = thread * nH_;
        for (unsigned int i = 0; i < nH_; i++)
          tempScd_[i] += tempScd_[offset+i];
      }
    }
#   else
    void AllocateTempSpace()              { tempScd_.assign( nH_, 0 );  }
#   endif
    void IncrementTempBy(int i, double d) { tempScd_[i] += d;           }
    double TempSCD(int i)           const { return tempScd_[i];         }
#   ifdef MPI
    double* Sptr()       { return &sum_[0];    }
    double* S2ptr()      { return &sum2_[0];   }
    unsigned int* Nptr() { return &nvals_;     }
#   endif /* MPI */
  private:
    NameType name_;      ///< Carbon name
    Darray tempScd_;     ///< temporary space for calculating SCD values each frame.
    double sum_[3];      ///< Hold order param sum for each C-HX
    double sum2_[3];     ///< Hold order param sum^2 for each C-HX
    unsigned int nvals_; ///< # times this list has been updated, used to calc avg. from sum
    unsigned int nH_;    ///< Number of hydrogens
    bool init_;          ///< False if name has not yet been set.
};

/// Hold position information for a single lipid carbon site.
class Action_LipidOrder::CarbonSite {
  public:
    CarbonSite();
    /// Create site with carbon atom index, chain index, and position in chain
    CarbonSite(int,int,int);
    bool operator<(CarbonSite const& rhs) const { return (c_idx_ < rhs.c_idx_); }
    unsigned int NumH() const { return nH_; }
    int Cidx()          const { return c_idx_;  }
    int ChainIdx()      const { return chainIdx_; }
    int Position()      const { return position_; }
    int Hidx(int i)     const { return h_idx_[i]; }
    void AddHindex(int);
  private:
    int c_idx_;       ///< Carbon atom index.
    int chainIdx_;    ///< Chain type index
    int position_;    ///< Position in chain.
    int h_idx_[3];    ///< Hydrogen atom indices.
    unsigned int nH_; ///< Number of hydrogen atom indices.
};
#endif
