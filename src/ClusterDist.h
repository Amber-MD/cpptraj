#ifndef INC_CLUSTERDIST_H
#define INC_CLUSTERDIST_H
#include "SymmetricRmsdCalc.h"
#include "DataSet_Coords.h"
#include "DataSet_1D.h"
/// Abstract Base Class for Cluster centroid.
/** This class is a container for the cluster centroid type appropriate for
  * the data being clustered. For COORDS DataSets this is a frame, for other
  * DataSets this is just a number. Centroid classes must implement a Copy()
  * function.
  */
class Centroid { 
  public:
    virtual ~Centroid() {}
    virtual Centroid* Copy() = 0;
  // TODO: Should centroids remember number of frames that went into them?
  //       This would make it so FrameOpCentroid wouldnt require extra arg.
};
/// Cluster centroid for generic DataSet.
class Centroid_Num : public Centroid {
  public:
    Centroid_Num()           : cval_(0.0), sumx_(0.0), sumy_(0.0) {}
    Centroid_Num(double val, double x, double y) : cval_(val), sumx_(x), sumy_(y) {}
    Centroid* Copy() { return (Centroid*)new Centroid_Num(cval_, sumx_, sumy_); }
    friend class ClusterDist_Num;
  private:
    double cval_;
    double sumx_; // For storing periodic average
    double sumy_; // For storing periodic average
};
/// Cluster centroid for multiple DataSets
class Centroid_Multi : public Centroid {
  public:
    typedef std::vector<double> Darray;
    Centroid_Multi() {}
    Centroid_Multi(Darray const& val, Darray const& x, Darray const& y) :
      cvals_(val), Sumx_(x), Sumy_(y) {}
    Centroid* Copy() { return (Centroid*)new Centroid_Multi(cvals_, Sumx_, Sumy_); }
    friend class ClusterDist_Euclid;
  private:
    Darray cvals_;
    Darray Sumx_; // For storing periodic average
    Darray Sumy_; // For storing periodic average
};
/// Cluster Centroid for Coords DataSet.
class Centroid_Coord : public Centroid {
  public:
    Centroid_Coord() {}
    Centroid_Coord(Frame const& frame) : cframe_(frame) {}
    Centroid_Coord(int natom) : cframe_(natom) {}
    Centroid* Copy() { return (Centroid*)new Centroid_Coord(cframe_); }
    friend class ClusterDist_DME;
    friend class ClusterDist_RMS;
    friend class ClusterDist_SRMSD;
  private:
    Frame cframe_;
};
// -----------------------------------------------------------------------------
/// Abstract Base Class for Cluster distance calc.
class ClusterDist {
  public:
    enum CentOpType { ADDFRAME=0, SUBTRACTFRAME };
    /// Used to pass in absolute frame numbers for centroid calculations.
    typedef std::vector<int> Cframes;
    typedef Cframes::const_iterator Cframes_it;
    typedef std::vector<DataSet*> DsArray;
    virtual ~ClusterDist() {}
    /// \return distance between given frames.
    virtual double FrameDist(int, int) = 0;
    /// \return distance between given centroids.
    virtual double CentroidDist( Centroid*, Centroid* ) = 0;
    /// \return distance between given frame and centroid.
    virtual double FrameCentroidDist(int, Centroid* ) = 0;
    /// Calculate centroid from given frames.
    virtual void CalculateCentroid(Centroid*, Cframes const&) = 0;
    /// \return new centroid from given frames.
    virtual Centroid* NewCentroid(Cframes const&) = 0;
    /// \return copy of this ClusterDist
    virtual ClusterDist* Copy() = 0;
    /// Update centroid by performing given operation between given frame and centroid.
    virtual void FrameOpCentroid(int, Centroid*, double, CentOpType) = 0;
    /// \return string containing description of the distance metric
    virtual std::string Description() const = 0;
  protected:
    typedef double (*DistCalc)(double,double);
};
/// Cluster distance calc for generic DataSet
class ClusterDist_Num : public ClusterDist {
  public:
    ClusterDist_Num() : data_(0), dcalc_(0) {}
    ClusterDist_Num(DataSet*);
    double FrameDist(int, int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
    void FrameOpCentroid(int, Centroid*, double, CentOpType);
    ClusterDist* Copy() { return new ClusterDist_Num( *this ); }
    std::string Description() const;
  private:
    DataSet_1D* data_;
    DistCalc dcalc_;
};
/// Cluster distance calc using Euclid distance
class ClusterDist_Euclid : public ClusterDist {
  public:
    ClusterDist_Euclid() {}
    ClusterDist_Euclid(DsArray const&);
    double FrameDist(int, int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
    void FrameOpCentroid(int, Centroid*, double, CentOpType);
    ClusterDist* Copy() { return new ClusterDist_Euclid( *this ); }
    std::string Description() const;
  private:
    typedef std::vector<DataSet_1D*> D1Array;
    D1Array dsets_;
    typedef std::vector<DistCalc> DcArray;
    DcArray dcalcs_;
};
/// DME cluster distance calc for Coords DataSet.
class ClusterDist_DME: public ClusterDist {
  public:
    ClusterDist_DME() : coords_(0) {}
    ClusterDist_DME(DataSet*,AtomMask const&);
    double FrameDist(int, int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
    void FrameOpCentroid(int, Centroid*, double, CentOpType);
    ClusterDist* Copy() { return new ClusterDist_DME( *this ); }
    std::string Description() const;
  private:
    DataSet_Coords* coords_;
    AtomMask mask_;
    Frame frm1_;             ///< Temporary storage for frames from coords
    Frame frm2_;             ///< Temporary storage for frames from coords
};
/// RMS cluster distance calc for Coords DataSet.
class ClusterDist_RMS : public ClusterDist {
  public:
    ClusterDist_RMS() : coords_(0), nofit_(false), useMass_(false) {}
    ClusterDist_RMS(DataSet*,AtomMask const&,bool,bool);
    double FrameDist(int, int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    void FrameOpCentroid(int, Centroid*, double, CentOpType);
    Centroid* NewCentroid(Cframes const&);
    ClusterDist* Copy() { return new ClusterDist_RMS( *this ); }
    std::string Description() const;
  private:
    DataSet_Coords* coords_;
    AtomMask mask_;
    bool nofit_;
    bool useMass_;
    Frame frm1_;
    Frame frm2_;
};
/// Symmetry-corrected RMS distance calc for Coords DataSet.
class ClusterDist_SRMSD : public ClusterDist {
  public:
    ClusterDist_SRMSD() {}
    ClusterDist_SRMSD(DataSet*,AtomMask const&,bool,bool,int);
    double FrameDist(int, int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
    void FrameOpCentroid(int, Centroid*, double, CentOpType);
    ClusterDist* Copy() { return new ClusterDist_SRMSD( * this ); }
    std::string Description() const;
  private:
    DataSet_Coords* coords_;
    AtomMask mask_;
    SymmetricRmsdCalc SRMSD_;
    Frame frm1_;
    Frame frm2_;
}; 
#endif
