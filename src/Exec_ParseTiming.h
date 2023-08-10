#ifndef INC_EXEC_PARSETIMING_H
#define INC_EXEC_PARSETIMING_H
#include "Exec.h"
/// For extracting and parsing timing data from cpptraj output 
class Exec_ParseTiming : public Exec {
  public:
    Exec_ParseTiming();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ParseTiming(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<double> Darray;
    class RunTiming;
    typedef std::vector<RunTiming> RunArray;

    // Which variables to plot
    enum Xtype { X_INDEX=0, X_CORES };
    enum Ytype { Y_T_TOTAL=0, Y_T_TRAJPROC, Y_T_TRAJREAD, Y_T_ACTFRAME };
    // Groupo types. NOTE: Update GroupTypeStr in Exec_ParseTiming.cpp
    enum GroupType { GROUPBY_PREFIX = 0, GROUPBY_NAME, GROUPBY_KIND };

    static inline double YVAL(Ytype, RunTiming const&);

    RunTiming read_cpptraj_output(std::string const&);
    int create_output_set(RunArray const&, DataSetList&, DataFile*,
                          std::string const&, Dimension const&,
                          Xtype, Ytype) const;
    void write_to_file(CpptrajFile&, RunArray const&, Xtype, Ytype, double, double) const;
};
#endif
