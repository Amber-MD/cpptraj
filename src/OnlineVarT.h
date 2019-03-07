// -*- mode: c++; -*-

#ifndef INC_ONLINEVART_H
#define INC_ONLINEVART_H

#include <map>

/** Class that allows numerically stable calculation of mean/variance
  * via Welford's Online algorithm.
  */
template <class Float>
class Stats {
public:
  /// Default CONSTRUCTOR
  Stats() : n_(0.0), mean_(0.0), M2_(0.0) {}
# ifdef MPI
  /// CONSTRUCTOR taking N, mean, and M2
  Stats(Float n, Float mean, Float m2) : n_(n), mean_(mean), M2_(m2) {}
  Float M2()    const { return M2_; }; // Needed for MPI send 
# endif

  void accumulate(const Float x)
  {
    Float delta;

    n_++;
    delta = x - mean_;
    mean_ += delta / n_;
    M2_ += delta * (x - mean_);
  }

  Float mean() const { return mean_; };
  Float variance() const { 
    if (n_ < 2) return 0.0;
    return M2_ / (n_ - 1.0); 
  };
  Float nData() const { return n_; };
  /// Combine two averages and variances into this Stats
  void Combine(Stats<Float> const& rhs) {
    Float nX = n_ + rhs.n_;
    Float delta = rhs.mean_ - mean_;
    // Combined mean
    mean_ = mean_ + delta * (rhs.n_ / nX);
    // Combined variance
    M2_ = M2_ + rhs.M2_ + ((delta*delta) * ((n_*rhs.n_) / nX));
    n_ = nX; 
  }
//# ifdef MPI
  
  //Float M2()    const { return M2_; }; // Needed for MPI reduce
  //void SetVals(double m0, double m2, double n) {
  //  n_ = n; mean_ = m0; M2_ = m2;
  //}
//# endif
private:
  Float n_;
  Float mean_;
  Float M2_;
};


/** Class that allows numerically stable accumulation of bin values
  * (average and variance) using Welford's Online algorithm. Intended
  * for use in e.g. histograms. Both Key and Value are of primitive type
  * \param Key should be of integer type
  *\param Value should be of floating point type
  */
template <typename Key, typename Value>
class StatsMap {
public:
  StatsMap() :
    n_(0.0), min_(0), max_(0)
  {}

  typedef typename std::map<Key,Value>::iterator iterator;
  typedef typename std::map<Key,Value>::const_iterator const_iterator;

  void accumulate(std::map<Key,Value> a)
  {
    Value delta;
    Key min, max;


    // FIXME: This is ugly but we must make sure that _all_ indices from min to
    //        max are generated so that the mean can be calculated for every
    //        iteration.

    min = a.begin()->first;
    max = a.rbegin()->first;

    if (min < min_)
      min_ = min;
    
    if (max > max_)
      max_ = max;

    n_++;

    for (Key i = min_; i <= max_; i++) {
      delta = a[i] - mean_[i];
      mean_[i] += delta / n_;
      M2_[i] += delta * (a[i] - mean_[i]);
    }
  }

  Value mean(Key i) { return mean_[i]; };
  Value variance(Key i) {
    if (n_ < 2) return 0.0;
    return M2_[i] / (n_ - 1.0); 
  };

  iterator mean_begin()             { return mean_.begin(); }
  iterator mean_end()               { return mean_.end();   }
  const_iterator mean_begin() const { return mean_.begin(); }
  const_iterator mean_end()   const { return mean_.end();   }


  // not really the variance so will have to be divided by n - 1
  iterator variance_begin()             { return M2_.begin(); }
  iterator variance_end()               { return M2_.end();   }
  const_iterator variance_begin() const { return M2_.begin(); }
  const_iterator variance_end()   const { return M2_.end();   }

  Value nData() const { return n_; };

  bool empty() const { return n_ < 1; }

private:
  Value n_;
  Key min_, max_;
  std::map<Key,Value> mean_;
  std::map<Key,Value> M2_;
};
#endif
