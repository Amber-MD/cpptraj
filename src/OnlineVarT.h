// -*- mode: c++; -*-

#ifndef INC_ONLINEVART_H
#define INC_ONLINEVART_H

#include <map>


template <class Float>
class Stats {
public:
  Stats() :
    n_(0.0),
    mean_(0.0),
    M2_(0.0)
  {}

  void accumulate(const Float x)
  {
    Float delta;

    n_++;
    delta = x - mean_;
    mean_ += delta / n_;
    M2_ += delta * (x - mean_);
  }

  Float mean() const { return mean_; };
  Float variance() const { return M2_ / (n_ - 1); };
  Float nData() const { return n_; };

private:
  Float n_;
  Float mean_;
  Float M2_;
};


// statistic map: both Key and Value are of primitive type
// Key should be of integer type
// Value should be of floating point type

template <typename Key, typename Value>
class StatsMap {
public:
  StatsMap() :
    n_(0.0)
  {}

  typedef typename std::map<Key,Value>::iterator iterator;

  void accumulate(std::map<Key,Value> a)
  {
    Key i;
    Value delta;

    n_++;

    for (iterator it = a.begin(); it != a.end(); it++) {
      i = it->first;
      delta = a[i] - mean_[i];
      mean_[i] += delta / n_;
      M2_[i] += delta * (a[i] - mean_[i]);
    }
  }

  Value mean(Key i) { return mean_[i]; };
  Value variance(Key i) { return M2_[i] / (n_ - 1); };

  iterator mean_begin() { return mean_.begin(); };
  iterator mean_end() { return mean_.end(); };

  // not really the variance so will have to be divided by n - 1
  iterator variance_begin() { return M2_.begin(); };
  iterator variance_end() { return M2_.end(); };

  Value nData() const { return n_; };

private:
  Value n_;
  std::map<Key,Value> mean_;
  std::map<Key,Value> M2_;
};
#endif
