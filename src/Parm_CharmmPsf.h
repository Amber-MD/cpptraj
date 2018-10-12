#ifndef INC_PARM_CHARMMPSF_H
#define INC_PARM_CHARMMPSF_H
#include "ParmIO.h"
#include "ParameterSet.h"
class Parm_CharmmPsf : public ParmIO {
  public :
    Parm_CharmmPsf();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_CharmmPsf(); }
    bool ID_ParmFormat(CpptrajFile&);
    static void ReadHelp();
    int processReadArgs(ArgList&);
    int ReadParm(FileName const&, Topology&);
    int WriteParm(FileName const&, Topology const&);
    static void WriteHelp();
    int processWriteArgs(ArgList&);
  private:
    static inline int FindTag(char*, const char*, CpptrajFile&);
    int ReadDihedrals(CpptrajFile&, int, const char*, Topology&) const;

    inline void WriteSectionHeader(CpptrajFile&, const char*, int) const;

    ParameterSet params_;
    bool extfmt_; ///< (write) If true use extended format
    bool cheq_;   ///< (write) If true include extra columns for polarization (CHarge EQuilibration)
    bool xplor_;  ///< (write) If true use XPLOR format PSF
};
#endif
