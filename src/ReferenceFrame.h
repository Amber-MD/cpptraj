#ifndef INC_REFERENCEFRAME_H
#define INC_REFERENCEFRAME_H
#include "DataSet_Coords_REF.h"
/// Wrapper around DataSet_Coords_REF DataSet.
/** Intended as a non-modifiable holder for the DataSet_Coords_REF DataSet
  * that also can hold an error status.
  */
class ReferenceFrame {
  public:
    ReferenceFrame() : ref_(0), err_(0) {}
    ReferenceFrame(int err) : ref_(0), err_(err) {}
    ReferenceFrame(DataSet_Coords_REF* ds) : ref_(ds), err_(0) {}
    ReferenceFrame(const ReferenceFrame& rhs) : ref_(rhs.ref_), err_(rhs.err_) {}
    ReferenceFrame& operator=(ReferenceFrame const& rhs) {
      if (this == &rhs) return *this;
      ref_ = rhs.ref_;
      err_ = rhs.err_;
      return *this;
    }
    Frame const& Coord()        const { return ref_->RefFrame();  }
    Topology const& Parm()      const { return ref_->Top();       }
    Topology* ParmPtr()               { return ref_->TopPtr();    } // FIXME deprecate this
    bool error()                const { return err_ != 0;         }
    bool empty()                const { return ref_ == 0;         }
    DataSet_Coords_REF* RefPtr()const { return ref_;              }
    CoordinateInfo const& CoordsInfo() const { return ref_->CoordsInfo(); }
    /// \return base file name, or if that is empty data set name.
    std::string const& RefName() const{
      if (ref_->Meta().Fname().empty())
        return ref_->Meta().Name();
      else
        return ref_->Meta().Fname().Base();
    }
    const char* refName() const { return RefName().c_str() ; }
  private:
    DataSet_Coords_REF* ref_; ///< Reference coords DataSet from e.g. DataSetList.
    int err_;                 ///< Error status.
};
#endif
