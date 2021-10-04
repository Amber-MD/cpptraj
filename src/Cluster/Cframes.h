#ifndef INC_CLUSTER_CFRAMES_H
#define INC_CLUSTER_CFRAMES_H
#include <vector>
#include <cstddef> // size_t
namespace Cpptraj {
namespace Cluster {

/// Used to hold frames to cluster (may not be all frames).
class Cframes {
    typedef std::vector<int> Iarray;
  public:
    Cframes() {}
    Cframes(std::size_t n, int i) : frames_(n, i) {}

    std::size_t size()               const { return frames_.size(); }
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
    void assign(std::size_t n, int i)      { frames_.assign(n, i);   }
    void clear()                           { frames_.clear();        }
    void reserve(std::size_t n)            { frames_.reserve(n);     }
    Iarray const& Data()             const { return frames_;         }
    bool empty()                     const { return frames_.empty(); }

    bool HasFrame(int) const;
    int FrameIdx(int) const;
    void Insert(Cframes const& rhs) { frames_.insert(frames_.end(), rhs.begin(), rhs.end()); }
    void Remove(int);
    void Sort();
  private:
    Iarray frames_;    ///< Frames to cluster.
};

typedef Cframes::const_iterator Cframes_it;

}
}
#endif
