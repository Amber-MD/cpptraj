#ifndef INC_STRUCTURE_RESSTATARRAY_H
#define INC_STRUCTURE_RESSTATARRAY_H
namespace Cpptraj {
namespace Structure {
/// Used to keep track of residue status during build
class ResStatArray {
  public:
    enum Type { UNKNOWN = 0,
                       VALIDATED,
                       SUGAR_UNRECOGNIZED_LINK_RES,
                       SUGAR_UNRECOGNIZED_LINKAGE,
                       SUGAR_NO_LINKAGE,
                       SUGAR_NO_CHAIN_FOR_LINK,
                       SUGAR_NAME_MISMATCH,
                       //SUGAR_MISSING_C1X,
                       SUGAR_SETUP_FAILED };

    /// CONSTRUCTOR
    ResStatArray() {}

    Type const& operator[](int idx) const { return resStat_[idx]; }

    Type& operator[](int idx) { return resStat_[idx]; }
  private:
    typedef std::vector<Type> Rarray;
    Rarray resStat_; ///< Status of each residue in a topology
};
}
}
#endif
