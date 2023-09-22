#ifndef INC_EXEC_CRDTRANSFORM_H
#define INC_EXEC_CRDTRANSFORM_H
#include "Exec.h"
#include "ExtendedSimilarity.h"
/// Used to transform/condition coordinates 
class Exec_CrdTransform : public Exec {
  public:
    Exec_CrdTransform() : Exec(COORDS), debug_(0) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CrdTransform(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    enum CriterionType { COMP_SIM = 0,   ///< Remove most dissimilar objects based on complement similarity
                         SIM_TO_MEDOID,  ///< Remove most dissimilar objects based on similarity to medoid.
                         NO_CRITERION };

    static const char* CriterionStr_[];

    int iterativeRmsRefinement(AtomMask const&, bool, double,
                               DataSet_Coords*, DataSet_Coords*) const;
    int normalizeCoords(DataSet_Coords*, DataSet_Coords*) const;
    int trimOutliers(int, double, ExtendedSimilarity::MetricType, CriterionType,
                     DataSet_Coords*, DataSet_Coords*) const;

    int debug_;
};
#endif
