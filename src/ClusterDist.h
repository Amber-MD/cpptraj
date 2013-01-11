#ifndef INC_CLUSTERDIST_H
#define INC_CLUSTERDIST_H
#include <list>
#include "Frame.h"
#include "DataSet_Coords.h"
#include "ClusterMatrix.h"
/// Abstract Base Class for Cluster centroid.
class Centroid { 
  public:
    virtual ~Centroid() {}
    virtual Centroid* Copy() = 0;
};
/// Cluster centroid for generic DataSet.
class Centroid_Num : public Centroid {
  public:
    Centroid_Num()           : cval_(0.0) {}
    Centroid_Num(double val) : cval_(val) {}
    Centroid* Copy() { return (Centroid*)new Centroid_Num(cval_); }
    friend class ClusterDist_Num;
  private:
    double cval_;
};
/// Cluster centroid for mulitple DataSets
class Centroid_Multi : public Centroid {
  public:
    Centroid_Multi() {}
    Centroid_Multi(std::vector<double> const& val) : cvals_(val) {}
    Centroid* Copy() { return (Centroid*)new Centroid_Multi(cvals_); }
    friend class ClusterDist_Euclid;
  private:
    std::vector<double> cvals_;
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
  private:
    Frame cframe_;
};

/// Abstract Base Class for Cluster distance calc.
/** The pairwise-distance calculation is here to make COORDS DataSet calcs 
  * more efficient; otherwise they would have to copy frame1 coords each time 
  * as well as always track memory for frame2.
  */
class ClusterDist {
  public:
    typedef std::list<int> Cframes;
    typedef Cframes::const_iterator Cframes_it;
    typedef std::vector<DataSet*> DsArray;
    virtual ~ClusterDist() {}
    virtual ClusterMatrix PairwiseDist(int) = 0;
    virtual double CentroidDist( Centroid*, Centroid* ) = 0;
    virtual double FrameCentroidDist(int, Centroid* ) = 0;
    virtual void CalculateCentroid(Centroid*, Cframes const&) = 0;
    virtual Centroid* NewCentroid(Cframes const&) = 0;
  protected:
    typedef double (*DistCalc)(double,double);
};
/// Cluster distance calc for generic DataSet
class ClusterDist_Num : public ClusterDist {
  public:
    ClusterDist_Num() : data_(0), dcalc_(0) {}
    ClusterDist_Num(DataSet*);
    ClusterMatrix PairwiseDist(int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
  private:
    DataSet* data_;
    DistCalc dcalc_;
};
/// Cluster distance calc using Euclid distance
class ClusterDist_Euclid : public ClusterDist {
  public:
    ClusterDist_Euclid() {}
    ClusterDist_Euclid(DsArray const&);
    ClusterMatrix PairwiseDist(int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
  private:
    DsArray dsets_;
    typedef std::vector<DistCalc> DcArray;
    DcArray dcalcs_;
};
/// DME cluster distance calc for Coords DataSet.
class ClusterDist_DME: public ClusterDist {
  public:
    ClusterDist_DME() : coords_(0) {}
    ClusterDist_DME(DataSet*,std::string const&);
    ClusterMatrix PairwiseDist(int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
  private:
    DataSet_Coords* coords_;
    AtomMask mask_;
    Frame frm1_;             ///< Temporary storage for frames from coords
};
/// RMS cluster distance calc for Coords DataSet.
class ClusterDist_RMS : public ClusterDist {
  public:
    ClusterDist_RMS() : coords_(0), nofit_(false), useMass_(false) {}
    ClusterDist_RMS(DataSet*,std::string const&,bool,bool);
    ClusterMatrix PairwiseDist(int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
  private:
    DataSet_Coords* coords_;
    AtomMask mask_;
    bool nofit_;
    bool useMass_;
    Frame frm1_;
};
#endif
