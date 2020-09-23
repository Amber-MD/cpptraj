#ifndef INC_SEGMENT_H
#define INC_SEGMENT_H
//namespace Cpptraj {
/// Define a contiguous segment of atoms
class Segment {
  public:
    /// CONSTRUCTOR
    Segment() : begin_(-1), end_(-1) {}
    /// CONSTRUCTOR - Segment size of 1
    Segment(int idx) : begin_(idx), end_(idx+1) {}
    /// CONSTRUCTOR - Beginning and end
    Segment(int beg, int end) : begin_(beg), end_(end) {}

    int Begin() const { return begin_; }
    int End()   const { return end_; }
    int Size()  const { return (end_ - begin_); }

    /// Sort by beginning atom, then ending atom.
    bool operator<(const Segment& rhs) const {
      if (begin_ == rhs.begin_) {
        return (end_ < rhs.end_);
      } else {
        return (begin_ < rhs.begin_);
      }
    }

    /// Update the end of the segment
    void UpdateEnd(int idx) { end_ = idx + 1; }
  private:
    int begin_;
    int end_;
};
//}
#endif
