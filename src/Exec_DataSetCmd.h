#ifndef INC_EXEC_DATASETCMD_H
#define INC_EXEC_DATASETCMD_H
#include "Exec.h"
/// Process DataSet-specific command
class Exec_DataSetCmd : public Exec {
  public:
    Exec_DataSetCmd() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_DataSetCmd(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    enum CriterionType { UNKNOWN_C = 0, AVERAGE, SIZE, SMODE, STYPE, N_C };
    static const char* CriterionKeys[];
    enum SelectType { UNKNOWN_S = 0, EQUAL, NOT_EQUAL, LESS_THAN, GREATER_THAN, BETWEEN, OUTSIDE, N_S };
    struct SelectPairType {
      SelectType type_;
      const char* key_;
    };
    static SelectPairType SelectKeys[];

    static void Help_ModifyPoints();
    RetType ModifyPoints(CpptrajState&, ArgList&, bool);
    RetType VectorCoord(CpptrajState&, ArgList&);
    RetType ChangeOutputFormat(CpptrajState const&, ArgList&);
    RetType Remove(CpptrajState&, ArgList&);
    RetType MakeXY(CpptrajState&, ArgList&);
    RetType Make2D(CpptrajState&, ArgList&);
    RetType Filter(CpptrajState&, ArgList&);
    RetType Concatenate(CpptrajState&, ArgList&);
    static void Help_ChangeDim();
    RetType ChangeDim(CpptrajState const&, ArgList&);
    RetType ChangeModeType(CpptrajState const&, ArgList&);
    RetType InvertSets(CpptrajState&, ArgList&);
};
#endif
