#ifndef INC_ACTION_CLUSTERDIHEDRAL_H
#define INC_ACTION_CLUSTERDIHEDRAL_H
#include "Action.h"
class DCnode {
  public:
    DCnode() : count_(0) {}
    DCnode(const DCnode& rhs) : BinIDs_(rhs.BinIDs_), count_(rhs.count_) {}
    DCnode& operator=(const DCnode& rhs) {
      if (this==&rhs) return *this;
      BinIDs_ = rhs.BinIDs_;
      count_ = rhs.count_;
      return *this;
    }
    bool operator<(const DCnode& rhs) { return (count_ < rhs.count_); }
    bool operator>(const DCnode& rhs) { return (count_ > rhs.count_); }
    bool operator==(const DCnode& rhs) { return (count_ == rhs.count_); }
  private:
    std::vector<int> BinIDs_;
    int count_;
};

class DCmask {
  public: 
    DCmask() : a1_(0), a2_(0), a3_(0), a4_(0), bins_(0) {}
    DCmask(int a1, int a2, int a3, int a4, int bins) :
           a1_(a1), a2_(a2), a3_(a3), a4_(a4), bins_(bins) {}
  private:
    int a1_;
    int a2_;
    int a3_;
    int a4_;
    int bins_;
};

class Action_ClusterDihedral : public Action {
  public:
    Action_ClusterDihedral();
  private:
    int ReadDihedrals(std::string const&);    

    int init();

    std::vector<DCnode> dcarray_;  ///< Hold counts for each bin# combo.
    std::vector<int> clusterNums_; ///< Hold DC # for each frame.
    std::vector<DCmask> DCmasks_;  ///< Hold 4 atom mask for each dihedral
    int phibins_;
    int psibins_;
    int CUT_;
    std::string outfile_;
    std::string framefile_; // filenames[1]
    std::string infofile_;  // filenames[2]
    std::string cvtfile_;   // filenames[3]i
    AtomMask mask_;
};
#endif
