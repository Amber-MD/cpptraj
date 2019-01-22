#ifndef INC_CLUSTER_CFRAMES_H
#define INC_CLUSTER_CFRAMES_H
#include <vector>
namespace Cpptraj {
namespace Cluster {

/// Used to hold frames to cluster (may not be all frames).
class Cframes {
    typedef std::vector<int> Iarray;
  public:
    enum SieveType { NONE=0, REGULAR, RANDOM };

    Cframes()                : type_(NONE), sieve_(1) {}
    Cframes(size_t n, int i) : frames_(n, i), type_(NONE), sieve_(1) {}

    size_t size()                    const { return frames_.size(); }
    typedef Iarray::const_iterator const_iterator;
    const_iterator begin()           const { return frames_.begin(); }
    const_iterator end()             const { return frames_.end();   }
/*    typedef Iarray::iterator iterator;
    iterator begin()             { return frames_.begin(); }
    iteratir end() */
    int front()                      const { return frames_.front(); }
    int operator[](unsigned int idx) const { return frames_[idx];    }
    int& operator[](unsigned int idx)      { return frames_[idx];    }
    void push_back(int i)                  { frames_.push_back( i ); }
    void assign(size_t n, int i)           { frames_.assign(n, i);   }

    bool HasFrame(int) const;
    void Insert(Cframes const& rhs) { frames_.insert(frames_.end(), rhs.begin(), rhs.end()); }
    void Remove(int);
    void Sort();

    int SetFramesToCluster(int, size_t, int);
  private:
    void DetermineTypeFromSieve(int);

    Iarray frames_;
    SieveType type_;
    int sieve_;
};

typedef Cframes::const_iterator Cframes_it;

}
}
#endif
