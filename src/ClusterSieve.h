#ifndef INC_CLUSTERSIEVE_H
#define INC_CLUSTERSIEVE_H
#include <vector>
#include <cstddef>
/// Used to map actual frame numbers to ClusterMatrix internal indices.
class ClusterSieve {
  public:
    enum SieveType { NONE=0, REGULAR, RANDOM };
    typedef std::vector<int> SievedFrames;
    ClusterSieve();
    void Clear();
    /// Setup no sieve, regular sieve, or random sieve.
    int SetSieve(int, size_t, int);
    /// Setup sieve from array: 'T'=sieved, 'F'=not sieved.
    int SetSieve(int, std::vector<char> const&);
    /// \return array of frames that will actually be clustered.
    SievedFrames const& Frames()     const { return idxToFrame_; }
    /// \return size of data in bytes
    size_t DataSize() const;
    /// \return an array index corresponding to a sieved frame.
    inline int FrameToIdx(int frame) const { return frameToIdx_[frame]; }
    /// \return frame number correspnding to array index.
    inline int IdxToFrame(int idx)   const { return idxToFrame_[idx];   }
    /// \return Original max number of frames.
    inline size_t MaxFrames()        const { return frameToIdx_.size(); }
    /// \return Actual number of frames after sieving.
    inline int ActualNframes()       const { return actualNframes_;     }
    /// \return Sieve value.
    inline int Sieve()               const { return sieve_;             }
    /// \return Sieve type.
    inline SieveType Type()          const { return type_;              }
  private:
    inline void DetermineTypeFromSieve(int);
    void MakeIdxToFrame();
    SieveType type_;              ///< Sieve type.
    int sieve_;                   ///< Sieve value; > 1 is regular, < -1 is random.
    int actualNframes_;           ///< Actual number of frames after sieving.
    std::vector<int> frameToIdx_; ///< Frame number to matrix index; -1 if frame was sieved out.
    std::vector<int> idxToFrame_; ///< Matrix index to frame number (size: actualNframes_).
};
#endif
