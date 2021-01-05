#ifndef INC_ACTION_HYDROGENBOND_H
#define INC_ACTION_HYDROGENBOND_H
#include <map>
#include "Action.h"
#include "ImageOption.h"
#include "DataSet_integer.h"
#ifdef TIMER
# include "Timer.h"
#endif
class Action_HydrogenBond : public Action {
  public:
    Action_HydrogenBond();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_HydrogenBond(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
#   ifdef MPI
    int SyncAction();
    Parallel::Comm trajComm_;
#   endif
    void Print();

    typedef std::vector<int> Iarray;

    class Site;
    class Hbond;
    class Bridge;

    inline double Angle(const double*, const double*, const double*, Box const&) const;
    inline int UU_Set_Idx(int,int) const;
    inline DataSet_integer* UUset(int,int,int);
    void AddUU(double,double,int,int,int,int);
    void AddUV(double,double,int,int,int,int,bool);
    void CalcSiteHbonds(int,double,Site const&,const double*,int,const double*,
                        Frame const&, int&);
    void CalcSolvHbonds(int,double,Site const&,const double*,int,const double*,
                        Frame const&, int&, bool);
    /// Update all hydrogen bond time series
    void UpdateSeries();
    /// Determine memory usage from # hbonds and time series
    std::string MemoryUsage(size_t, size_t, size_t) const;
#   ifdef MPI
    static std::vector<int> GetRankNhbonds(int,Parallel::Comm const&);
#   endif

    typedef std::vector<Site> Sarray;
    typedef std::pair<int,int> Hpair;
    typedef std::map<Hpair,Hbond> UUmapType;
    typedef std::map<int,Hbond> UVmapType;
    typedef std::map< int,std::set<int> > RmapType;
    typedef std::map< std::set<int>,Bridge > BmapType;
    typedef std::vector<Hbond> Harray;
    typedef std::map<int,int> IdxMapType;

    Sarray Both_;         ///< Array of donor sites that can also be acceptors
    Iarray Acceptor_;     ///< Array of acceptor-only atom indices
    Sarray SolventSites_; ///< Array of solvent donor/acceptor sites

    UUmapType UU_Map_;        ///< Map solute donorH/acceptor pair to UU hbond
    UVmapType UV_Map_;        ///< Map solute donorH or solute acceptor to UV hbond
    RmapType solvent2solute_; ///< Map solvent res # to solute residues it is bound to each frame
    BmapType BridgeMap_; ///< Map residues involved in bridging to # frames bridge present
    IdxMapType DidxMap_; ///< Map solute hydrogen donor atom # to index (series only)
    IdxMapType AidxMap_; ///< Map solute acceptor atom # to index (series only).

#   ifdef _OPENMP
    std::vector<Harray> thread_HBs_; ///< Hold hbonds found by each thread each frame.
#   endif

    std::string hbsetname_; ///< DataSet name
    AtomMask DonorMask_;
    AtomMask DonorHmask_;
    AtomMask AcceptorMask_;
    AtomMask SolventDonorMask_;
    AtomMask SolventAcceptorMask_;
    AtomMask Mask_;
    ImageOption imageOpt_;       ///< Used to determine if imaging should be performed
#   ifdef TIMER
    Timer t_action_;
    Timer t_uu_;
    Timer t_uv_;
    Timer t_bridge_;
#   endif
    Topology* CurrentParm_;  ///< Used to set atom/residue labels
    DataSetList* masterDSL_; ///< Used to add series data
    DataSet* NumHbonds_;     ///< Hold # UU hbonds per frame.
    DataSet* NumSolvent_;    ///< Hold # UV hbonds per frame.
    DataSet* NumBridge_;     ///< Hold # solute-solvent bridges per frame.
    DataSet* BridgeID_;      ///< Hold info on each bridge per frame.
    DataFile* UUseriesout_;  ///< File to write UU time series to.
    DataFile* UVseriesout_;  ///< File to write UN time series to.
    CpptrajFile* avgout_;    ///< File to write UU averages to.
    CpptrajFile* solvout_;   ///< File to write UV averages to.
    CpptrajFile* bridgeout_; ///< File to write bridge totals to.
    double dcut2_;           ///< Heavy atom distance cutoff squared.
    double acut_;            ///< Angle cutoff in radians.
    unsigned int bothEnd_;   ///< Index in Both_ where donor-only sites begin
    int Nframes_;            ///< Number of frames action has been active
    int debug_;
    bool series_;             ///< If true track hbond time series.
    bool seriesUpdated_;      ///< If false hbond time series need to be finished.
    bool useAtomNum_;         ///< If true include atom numbers in labels/legends
    bool noIntramol_;         ///< If true ignore intramolecular hydrogen bonds/bridges.
    bool hasDonorMask_;       ///< If true a donor mask was specified.
    bool hasDonorHmask_;      ///< If true a donor H mask was specified.
    bool hasAcceptorMask_;    ///< If true an acceptor mask was specified.
    bool hasSolventDonor_;    ///< If true a solvent donor mask was specified.
    bool hasSolventAcceptor_; ///< If true a solvent acceptor mask was specified.
    bool calcSolvent_;        ///< If true solute-solvent hbonds and bridges will be calcd.
    bool bridgeByAtom_;       ///< If true determine bridging by atom.
    // TODO replace with class
    typedef std::pair< std::set<int>,int > Bpair;
    /// \return true if 1) p0 frames > p1 frames, 2) p0 res < p1 res
    struct bridge_cmp {
      inline bool operator()(Bpair const& p0, Bpair const& p1) const
      {
        if (p0.second == p1.second)
          return (p0.first < p1.first);
        else
          return (p0.second > p1.second);
      }
    };
};

// ----- CLASSES ---------------------------------------------------------------
/// Potential hydrogen bond site. Can be either donor or donor/acceptor.
class Action_HydrogenBond::Site {
  public:
    Site() : idx_(-1) {}
    /// Solute site - heavy atom, hydrogen atom
    Site(int d, int h) : hlist_(1,h), idx_(d) {}
    /// Solute site - heavy atom, list of hydrogen atoms
    Site(int d, Iarray const& H) : hlist_(H), idx_(d) {}
    /// \return heavy atom index
    int Idx() const { return idx_; }
    /// \return number of hydrogen indices
    unsigned int n_hydrogens()      const { return hlist_.size(); }
    /// \return true if site is an ion (D atom == H atom)
    bool IsIon() const { return (hlist_.size()==1 && hlist_[0] == idx_); }
    /// \return iterator to beginning of hydrogen indices
    Iarray::const_iterator Hbegin() const { return hlist_.begin(); }
    /// \return iterator to end of hydrogen indices
    Iarray::const_iterator Hend()   const { return hlist_.end(); }
  private:
    Iarray hlist_; ///< List of hydrogen indices
    int idx_;      ///< Heavy atom index
};

/// Track specific hydrogen bond.
class Action_HydrogenBond::Hbond {
  public:
    Hbond() : dist_(0.0), angle_(0.0), data_(0), A_(-1), H_(-1), D_(-1), frames_(0) {}
    /// New hydrogen bond
    Hbond(double d, double a, DataSet_integer* s, int ia, int ih, int id) :
      dist_(d), angle_(a), data_(s), A_(ia), H_(ih), D_(id), frames_(1) {}
#   ifdef _OPENMP
    /// Just record that hbond exists
    Hbond(double d, double a, int ia, int ih, int id) :
      dist_(d), angle_(a), data_(0), A_(ia), H_(ih), D_(id), frames_(0) {}
    /// This version is for UV hbonds; a 1 in frames_ indicates soluteDonor
    Hbond(double d, double a, int ia, int ih, int id, int sd) :
      dist_(d), angle_(a), data_(0), A_(ia), H_(ih), D_(id), frames_(sd) {}
#   endif
    double Dist()  const { return dist_;   }
    double Angle() const { return angle_;  }
    int Frames()   const { return frames_; }
    int A()        const { return A_;      }
    int H()        const { return H_;      }
    int D()        const { return D_;      }
#   ifdef MPI
    DataSet_integer* Data() const { return data_; }
    /// New hydrogen bond with given # frames
    Hbond(double d, double a, DataSet_integer* s, int ia, int ih, int id, int n) :
      dist_(d), angle_(a), data_(s), A_(ia), H_(ih), D_(id), frames_(n) {}
    /// Update distance/angle/number frames
    void Combine(double d, double a, int n) {
      dist_ += d;
      angle_ += a;
      frames_ += n;
    }
#   endif
    /// First sort by frames (descending), then distance (ascending).
    bool operator<(const Hbond& rhs) const {
      if (frames_ == rhs.frames_)
        return (dist_ < rhs.dist_);
      else
        return (frames_ > rhs.frames_);
    }
    /// Update distance/angle/time series
    void Update(double d, double a, int f) {
      dist_ += d;
      angle_ += a;
      ++frames_;
      if (data_ != 0) data_->AddVal(f, 1);
    }
    void CalcAvg();
    void FinishSeries(unsigned int);
  private:
    static const int ZERO;
    double dist_;  ///< Used to calculate average distance of hydrogen bond
    double angle_; ///< Used to calculate average angle of hydrogen bond
    DataSet_integer* data_; ///< Hold time series data
    int A_; ///< Acceptor atom index
    int H_; ///< Hydrogen atom index
    int D_; ///< Donor atom index
    int frames_; ///< # frames this hydrogen bond has been present
};

/// Track solvent bridge between 2 or more solute residues.
class Action_HydrogenBond::Bridge {
  public:
    /// Constructor - new bridge
    Bridge() : data_(0), frames_(1) {}
#   ifdef MPI
    /// Constructor - new bridge with given # frames
    Bridge(int f) : data_(0), frames_(f) {}
    /// Increment number of frames
    void Combine(int n) { frames_ += n; }
#   endif
    int Frames() const { return frames_; }
    /// Update frames/time series
    void Update(int f) {
      ++frames_;
      if (data_ != 0) data_->AddVal(f, 1);
    }
    /// \return true if bridge has more frames than rhs.
    bool operator()(Bridge const& rhs) const {
      if (frames_ > rhs.frames_)
        return true;
      else if (frames_ < rhs.frames_)
        return false;
      else
        return false;
    }
  private:
    DataSet_integer* data_; ///< Hold time series data
    int frames_; ///< # frames this bridge has been present.
};
#endif
