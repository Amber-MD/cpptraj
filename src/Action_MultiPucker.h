#ifndef INC_ACTION_MULTIPUCKER_H
#define INC_ACTION_MULTIPUCKER_H
#include "Action.h"
#include "Pucker.h"
#include "Pucker_PuckerSearch.h"
#include "Range.h"
/// Automatically detect and calculate puckers within a residue range. 
class Action_MultiPucker : public Action {
  public:
    Action_MultiPucker();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_MultiPucker(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    static const double PERIOD_;                 ///< Pucker period in degrees (360)

    Cpptraj::Pucker::PuckerSearch puckerSearch_; ///< Used to search for puckers
    std::vector<DataSet*> data_;                 ///< Output DataSets, 1 per pucker
    std::vector<DataSet*> amp_;                  ///< Output amplitude DataSets, 1 per pucker
    std::vector<DataSet*> theta_;                ///< Output theta DataSets, 1 per pucker
    std::vector<Cpptraj::Pucker::Method> puckerMethods_; ///< Method to use for each pucker
    Range resRange_;                             ///< Residue range to search
    std::string dsetname_;                       ///< Output data set(s) name
    DataFile* outfile_;                          ///< File to write sets to
    DataFile* ampfile_;                          ///< File to write amplitude sets to
    DataFile* thetafile_;                        ///< File to write theta sets to
    DataSetList* masterDSL_;                     ///< Pointer to master DataSetList
    Cpptraj::Pucker::Method defaultMethod_;      ///< Which calculation method to use.
    double puckerMin_;                           ///< Min pucker value; set to 0 or -180
    double puckerMax_;                           ///< Max pucker value; set to 360 or 180
    double offset_;                              ///< Offset to add to pucker values.
    bool calc_amp_;                              ///< If true save amplitude as well
    bool calc_theta_;                            ///< If true save theta as well
};
#endif
