#ifndef INC_DATASET_COORDS_H
#define INC_DATASET_COORDS_H
#include "DataSet.h"
#include "Topology.h"
class DataSet_Coords : public DataSet {
  public:
    typedef std::vector<float> CRD;
    DataSet_Coords();
    ~DataSet_Coords();

    int Allocate(int);
    int Xmax() { return (int)(coords_.size() - 1); }
    int Size() { return (int)coords_.size(); }
    void Info();

    void AddFrame(Frame const&);
    int SetupTopMask(std::string const&, Topology&);
    int Natom()                    { return top_->Natom(); }
    const AtomMask& Mask()         { return mask_;         }
    const CRD& operator[](int idx) { return coords_[idx];  }
    Topology* Top()                { return top_;          }
  private:
    AtomMask mask_;                    ///< Mask
    Topology* top_;                    ///< Associated topology corresponding to mask.
    typedef std::vector<CRD> CRDarray;
    CRDarray coords_;                  ///< Array of coordinates.
};
#endif
