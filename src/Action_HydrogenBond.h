#ifndef INC_ACTION_HYDROGENBOND_H
#define INC_ACTION_HYDROGENBOND_H
#include <map>
#include "Action.h"
#include "ImageOption.h"
#include "DataSet_integer.h"
#include "OnlineVarT.h"
//#inc lude "CpptrajStdio.h" // DEBUG
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
    class bridgeSorter;
    /// Calcuate angle with optional imaging.
    inline double Angle(const double*, const double*, const double*, Box const&) const;
    /// \return solute-solute hydrogen bond index based on donor/acceptor atoms
    inline int UU_Set_Idx(int,int) const;
    /// \return solute-solute hydrogen bond time series DataSet.
    inline DataSet_integer* UUset(int,int,int);
    /// Add/update solute-solute hydrogen bond 
    void AddUU(double,double,int,int,int,int,int);
    /// Add/update solute-solvent hydrogen bond
    void AddUV(double,double,int,int,int,int,bool,int);
    /// Calculate UU hydrogen bonds between site and acceptor atom.
    void CalcSiteHbonds(int,double,Site const&,const double*,int,const double*,
                        Frame const&, int&, int);
    /// Calculate UV hydrogen bonds between solute/solvent site and solute/solvent acceptor atom.
    void CalcSolvHbonds(int,double,Site const&,const double*,int,const double*,
                        Frame const&, int&, bool, int);
    /// Update all hydrogen bond time series
    void UpdateSeries();
    /// Determine memory usage from # hbonds and time series
    std::string MemoryUsage(size_t, size_t, size_t) const;
    /// Print header for summary of parts
    static void summary_Parts_header(CpptrajFile*, unsigned int);
    /// Print summary of parts for hbond
    void summary_Parts(CpptrajFile*, Hbond const&) const;
#   ifdef MPI
    /// \return Array containing # hbonds on each rank
    static std::vector<int> GetRankNhbonds(int,Parallel::Comm const&);
    /// Flatten given Hbond into arrays
    static void HbondToArray(std::vector<double>&, std::vector<int>&, Hbond const&);
#   endif
    /// Ensure DataSet time series has at least N frames, fill out if necessary.
    static inline void FinishSeries(DataSet_integer*, unsigned int);

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
    Iarray splitFrames_;         ///< For calculating hydrogen bonds by parts
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
    DataFile* UVseriesout_;  ///< File to write UV time series to.
    DataFile* Bseriesout_;   ///< File to write bridge time series to.
    CpptrajFile* avgout_;    ///< File to write UU averages to.
    CpptrajFile* solvout_;   ///< File to write UV averages to.
    CpptrajFile* bridgeout_; ///< File to write bridge totals to.
    double dcut2_;           ///< Heavy atom distance cutoff squared.
    double acut_;            ///< Angle cutoff in radians.
    unsigned int bothEnd_;   ///< Index in Both_ where donor-only sites begin
    int Nframes_;            ///< Number of frames action has been active
    int debug_;
    bool series_;             ///< If true track hbond time series.
    bool Bseries_;            ///< If true track bridge time series.
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
// -----------------------------------------------------------------------------
/// Track specific hydrogen bond.
class Action_HydrogenBond::Hbond {
  public:
    Hbond() : dist_(0.0), angle_(0.0), data_(0), A_(-1), H_(-1), D_(-1), frames_(0) {}
    /// New hydrogen bond
    Hbond(DataSet_integer* s, int ia, int ih, int id, Iarray const& splits) :
      dist_(0), angle_(0), data_(s), A_(ia), H_(ih), D_(id), frames_(0)
    {
      if (!splits.empty()) {
        partsDist_.resize(splits.size()+1);
        partsAng_.resize(splits.size()+1);
      }
    }
#   ifdef _OPENMP
    /// Just record that hbond exists
    Hbond(double d, double a, int ia, int ih, int id) :
      dist_(d), angle_(a), data_(0), A_(ia), H_(ih), D_(id), frames_(0) {}
    /// This version is for UV hbonds; a 1 in frames_ indicates soluteDonor
    Hbond(double d, double a, int ia, int ih, int id, int sd) :
      dist_(d), angle_(a), data_(0), A_(ia), H_(ih), D_(id), frames_(sd) {}
#   endif
    double Dist()           const { return dist_;   }
    double Angle()          const { return angle_;  }
    int Frames()            const { return frames_; }
    int A()                 const { return A_;      }
    int H()                 const { return H_;      }
    int D()                 const { return D_;      }
    DataSet_integer* Data() const { return data_; }
    // Summary by parts
    unsigned int Nparts() const { return partsDist_.size(); }
    unsigned int PartFrames(unsigned int idx) const { return (unsigned int)partsDist_[idx].nData(); }
    double PartFrac(unsigned int idx, unsigned int Nframes) const { return partsDist_[idx].nData() / (double)Nframes; }
    Stats<double> const& PartDist(unsigned int idx)  const { return partsDist_[idx]; }
    Stats<double> const& PartAngle(unsigned int idx) const { return partsAng_[idx]; }
#   ifdef MPI
    /// CONSTRUCTOR - New hydrogen bond with given # frames
    Hbond(double d, double a, DataSet_integer* s, int ia, int ih, int id, int n) :
      dist_(d), angle_(a), data_(s), A_(ia), H_(ih), D_(id), frames_(n) {}
    /// Update distance/angle/number frames
    void Combine(double d, double a, int n) {
      dist_ += d;
      angle_ += a;
      frames_ += n;
    }
    /// Set up parts with given double array containing N, mean, and M2 for each part
    void SetupParts(unsigned int nparts, const double* dvals) {
      if (nparts == 0) return;
      partsDist_.clear();
      partsAng_.clear();
      partsDist_.reserve( nparts );
      partsAng_.reserve( nparts );
      const double* dptr = dvals;
      for (unsigned int idx = 0; idx != nparts; idx++, dptr += 6) {
        partsDist_.push_back( Stats<double>(dptr[0], dptr[1], dptr[2]) );
        partsAng_.push_back( Stats<double>(dptr[3], dptr[4], dptr[5]) );
      }
    }
    /// Update parts with given double array containing N, mean, and M2 for each part
    void CombineParts(unsigned int nparts, const double* dvals) {
      const double* dptr = dvals;
      for (unsigned int idx = 0; idx != nparts; idx++, dptr += 6) {
        partsDist_[idx].Combine( Stats<double>(dptr[0], dptr[1], dptr[2]) );
        partsAng_[idx].Combine( Stats<double>(dptr[3], dptr[4], dptr[5]) );
      }
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
    void Update(double distIn, double angIn, int fnum, Iarray const& splitFrames, int onum) {
      dist_ += distIn;
      angle_ += angIn;
      ++frames_;
      if (data_ != 0) data_->AddVal(fnum, 1);
      if (!splitFrames.empty()) {
        //mprintf("DEBUG: SPLIT: onum= %i partsDistSize= %zu splitFramesSize= %zu\n", onum, partsDist_.size(), splitFrames.size());
        // Find the correct part NOTE assumes onum never out of range
        unsigned int part = 0;
        while (part < splitFrames.size() && onum >= splitFrames[part]) part++;
        partsDist_[part].accumulate( distIn );
        partsAng_[part].accumulate( angIn );
      }
    }
    void CalcAvg();
  private:
    double dist_;  ///< Used to calculate average distance of hydrogen bond
    double angle_; ///< Used to calculate average angle of hydrogen bond
    DataSet_integer* data_; ///< Hold time series data
    int A_; ///< Acceptor atom index
    int H_; ///< Hydrogen atom index
    int D_; ///< Donor atom index
    int frames_; ///< # frames this hydrogen bond has been present
    std::vector<Stats<double>> partsDist_; ///< Hold avg. distance, split by parts
    std::vector<Stats<double>> partsAng_;  ///< Hold avg. angle, split by parts
};
// -----------------------------------------------------------------------------
/// Track solvent bridge between 2 or more solute residues.
class Action_HydrogenBond::Bridge {
  public:
    /// CONSTRUCTOR - new bridge
    Bridge(DataSet_integer* bds, Iarray const& splits) : data_(bds), frames_(0) {
      if (!splits.empty())
        partsFrames_.assign(splits.size()+1, 0);
    }
#   ifdef MPI
    /// Constructor - new bridge with given # frames
    Bridge(DataSet_integer* bds, int f) : data_(bds), frames_(f) {}
    /// Increment number of frames
    void Combine(int n) { frames_ += n; }
    /// Set up parts with given int array containing # frames for each part
    void SetupParts(unsigned int nparts, const int* ivals) {
      if (nparts == 0) return;
      partsFrames_.clear();
      partsFrames_.reserve( nparts );
      for (unsigned int idx = 0; idx != nparts; idx++)
        partsFrames_.push_back( ivals[idx] );
    }
    /// Update parts with given int array containing # frames for each part
    void CombineParts(unsigned int nparts, const int* ivals) {
      for (unsigned int idx = 0; idx != nparts; idx++)
        partsFrames_[idx] += ivals[idx];
    }
#   endif
    /// \return internal data set
    DataSet_integer* Data() const { return data_; }
    int Frames()            const { return frames_; }
    /// Update frames/time series
    void Update(int fnum, Iarray const& splitFrames, int onum) {
      ++frames_;
      if (data_ != 0) data_->AddVal(fnum, 1);
      if (!splitFrames.empty()) {
        // Find the correct part NOTE assumes onum never out of range
        unsigned int part = 0;
        while (part < splitFrames.size() && onum >= splitFrames[part]) part++;
        partsFrames_[part]++;
      }
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
    // Summary by parts
    unsigned int Nparts() const { return partsFrames_.size(); }
    unsigned int PartFrames(unsigned int idx) const { return (unsigned int)partsFrames_[idx]; }
  private:
    DataSet_integer* data_; ///< Hold time series data TODO
    int frames_;            ///< # frames this bridge has been present.
    Iarray partsFrames_;    ///< Hold # frames bridge present for each part.
};
#endif
