#ifndef INC_ACTION_HYDROGENBOND_H
#define INC_ACTION_HYDROGENBOND_H
#include <map>
#include "Action.h"
#include "ImagedAction.h"
#include "DataSet_integer.h"
#include "Timer.h"
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

    inline double Angle(const double*, const double*, const double*) const;
    void CalcSiteHbonds(int,double,Site const&,const double*,int,const double*,
                        Frame const&, int&);
    void CalcSolvHbonds(int,double,Site const&,const double*,int,const double*,
                        Frame const&, int&, bool);
    /// Update all hydrogen bond time series
    void UpdateSeries();
    /// Determine memory usage from # hbonds and time series
    std::string MemoryUsage(size_t, size_t) const;

    typedef std::vector<Site> Sarray;
    typedef std::pair<int,int> Hpair;
    typedef std::map<Hpair,Hbond> HBmapType;
    typedef std::map<int,Hbond> UVmapType;

    ImagedAction Image_; ///< Hold imaging info.
    Sarray Both_;     ///< Array of donor sites that can also be acceptors
    Iarray Acceptor_; ///< Array of acceptor-only atom indices
    Sarray SolventSites_;

    HBmapType UU_Map_;
    UVmapType UV_Map_;

    std::string hbsetname_;
    AtomMask DonorMask_;
    AtomMask DonorHmask_;
    AtomMask AcceptorMask_;
    AtomMask SolventDonorMask_;
    AtomMask SolventAcceptorMask_;
    AtomMask Mask_;
    Matrix_3x3 ucell_, recip_;
    Timer t_action_;
    Timer t_uu_;
    Timer t_uv_;
    Timer t_bridge_;
    Topology* CurrentParm_; ///< Used to set atom/residue labels
    DataSetList* masterDSL_;
    DataSet* NumHbonds_;
    DataSet* NumSolvent_;
    DataSet* NumBridge_;
    DataSet* BridgeID_;
    DataFile* UUseriesout_;
    DataFile* UVseriesout_;
    CpptrajFile* avgout_;
    CpptrajFile* solvout_;
    CpptrajFile* bridgeout_;
    double dcut2_;
    double acut_;
    unsigned int bothEnd_; ///< Index in Both_ where donor-only sites begin
    int Nframes_;        ///< Number of frames action has been active
    int debug_;
    bool series_;        ///< If true track hbond time series.
    bool seriesUpdated_; ///< If false hbond time series need to be finished.
    bool useAtomNum_;    ///< If true include atom numbers in labels/legends
    bool noIntramol_;
    bool hasDonorMask_;
    bool hasDonorHmask_;
    bool hasAcceptorMask_;
    bool hasSolventDonor_;
    bool calcSolvent_;
    bool hasSolventAcceptor_;
};

// ----- CLASSES ---------------------------------------------------------------
/// Potential hydrogen bond site. Can be either donor or donor/acceptor.
class Action_HydrogenBond::Site {
  public:
    Site() : idx_(-1), isV_(false) {}
    /// Solute site - heavy atom, hydrogen atom
    Site(int d, int h) : hlist_(1,h), idx_(d), isV_(false) {}
    /// Solute site - heavy atom, list of hydrogen atoms
    Site(int d, Iarray const& H) : hlist_(H), idx_(d), isV_(false) {}
    /// \return heavy atom index
    int Idx() const { return idx_; }
    /// \return iterator to beginning of hydrogen indices
    Iarray::const_iterator Hbegin() const { return hlist_.begin(); }
    /// \return iterator to end of hydrogen indices
    Iarray::const_iterator Hend()   const { return hlist_.end(); }
  private:
    Iarray hlist_; ///< List of hydrogen indices
    int idx_;      ///< Heavy atom index
    bool isV_;     ///< True if site is solvent
};

/// Track specific hydrogen bond.
class Action_HydrogenBond::Hbond {
  public:
    Hbond() : dist_(0.0), angle_(0.0), data_(0), A_(-1), H_(-1), D_(-1), frames_(0) {}
    /// New hydrogen bond
    Hbond(double d, double a, DataSet_integer* s, int ia, int ih, int id) :
      dist_(d), angle_(a), data_(s), A_(ia), H_(ih), D_(id), frames_(1) {}
    double Dist()  const { return dist_;   }
    double Angle() const { return angle_;  }
    int Frames()   const { return frames_; }
    int A()        const { return A_;      }
    int H()        const { return H_;      }
    int D()        const { return D_;      }
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
#endif
