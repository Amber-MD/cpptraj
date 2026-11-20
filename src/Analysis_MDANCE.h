#ifndef INC_ANALYSIS_MDANCE_H
#define INC_ANALYSIS_MDANCE_H
#include "Analysis.h"
#include "ExtendedSimilarity.h" // Has the same metrics
#include "TrajectoryFile.h"
#ifdef HAS_EIGEN
#include "Mdance/types.h"
#endif 
class DataSet_Coords;
/// Interface to MDANCE
/** MDANCE: Molecular Dynamics Analysis with N-aray Clustering Ensembles.
  */ 
class Analysis_MDANCE : public Analysis {
  public:
    Analysis_MDANCE();
    ~Analysis_MDANCE();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_MDANCE(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    static const char* kinitKeys_[];
    static const char* kinitStr_[];
#   ifdef HAS_EIGEN
    static const Cpptraj::Mdance::MD::KinitType kinitTypes_[];
#   endif
    typedef std::vector<int> Iarray;
    /// Hold frames/center for MDANCE clustering
    class MdCluster {
      public:
        MdCluster() : msd_(-1), rep_(-1) {}
        Iarray const& Frames() const { return frames_; }
        Frame const& Ctr() const { return ctr_; }
        unsigned int size() const { return frames_.size(); }
        int Rep() const { return rep_; }
        double MSD() const { return msd_; }

        void push_back(int i) { frames_.push_back( i ); }
        void SetCtr( Frame const& f ) { ctr_ = f; }
        void SetRep(int r) { rep_ = r; }
        void SetMSD(double m) { msd_ = m; }

        bool operator<(MdCluster const& rhs) const {
          if (frames_.size() > rhs.frames_.size())
            return true; // TODO absolute cluster number
          else
            return false;
        }
      private:
        Iarray frames_;
        Frame ctr_;
        double msd_;
        int rep_;
    };
    /// Hold MDANCE cluster results
    class ClusterArray {
      public:
        /// CONSTRUCTOR - take total # of things being clustered
        ClusterArray(unsigned int nframes) : nframes_(nframes) {}
        /// Reserve space for # of clusters
        void resize(unsigned int nclusters) { clusters_.resize(nclusters); }

        /// \return reference to specified cluster
        MdCluster& operator[](int cnum) { return clusters_[cnum]; }
        /// Set the DBI score
        void SetDBI(double d) { dbi_ = d; }
        /// Set the pseudo-F (Calinski-Harabasz) score
        void SetPSF(double p) { psf_ = p; }

        typedef std::vector<MdCluster>::iterator iterator;
        iterator begin() { return clusters_.begin(); }
        iterator end() { return clusters_.end(); }

        /// \return number of clusters
        unsigned int size() const { return clusters_.size(); }
        /// \return const reference to specified cluster
        MdCluster const& operator[](int cnum) const { return clusters_[cnum]; }
        /// \return DBI score
        double DBI() const { return dbi_; }
        /// \return pseudo-F (Calinkski-Harabasz) score
        double PSF() const { return psf_; }
        /// \return Number of frames clustered
        unsigned int Nframes() const { return nframes_; }

        typedef std::vector<MdCluster>::const_iterator const_iterator;
        const_iterator begin() const { return clusters_.begin(); }
        const_iterator end() const { return clusters_.end(); }

      private:
        std::vector<MdCluster> clusters_;
        double dbi_;
        double psf_;
        unsigned int nframes_;
    };

    void getClusterTrajArgs(ArgList&, const char*, const char*, std::string&,
                            TrajectoryFile::TrajFormatType&) const;
    void writeClusterTraj(ClusterArray const&) const;
    void writeCenterTraj(ClusterArray const&) const;
    void writeSummary(CpptrajFile&, ClusterArray const&) const;
    void writeInfo(CpptrajFile&, ClusterArray const&) const;
    void writeJson(CpptrajFile&, ClusterArray const&) const;

    int debug_;
    DataSet_Coords* coords_;
    int kClusters_;
    int kSeed_;
    int percentage_;
    int vthresh_;
    ExtendedSimilarity::MetricType metric_;
#   ifdef HAS_EIGEN
    Cpptraj::Mdance::MD::KinitType kinit_;
#   endif
    AtomMask mask_;

    DataSet* cnumvtime_;                        ///< Hold cluster number for each frame
    DataSet_Coords* centers_;                   ///< Hold cluster centers
    std::string clusterfile_;                   ///< Cluster trajectory base filename.
    std::string centerfile_;
    TrajectoryFile::TrajFormatType clusterfmt_; ///< Cluster trajectory format.
    TrajectoryFile::TrajFormatType centerfmt_;
    CpptrajFile* infofile_;
    CpptrajFile* summaryfile_;
    CpptrajFile* jsonfile_;
    bool sort_; ///< If true, sort by cluster population (default)
};
#endif
