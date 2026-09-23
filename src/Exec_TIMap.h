#ifndef INC_EXEC_TIMAP_H
#define INC_EXEC_TIMAP_H
#include "Exec.h"
#include "NameType.h"
#include "TIMatch.h"
#include <string>
#include <vector>
class Atom;
class DataSet;
class DataSetList;
class FileName;
class Frame;
class NameType;
class Topology;
/// Immediate command: write a residue library in template atom order; optional TI pair.
/** Command: timap (alias: timatch). See Help() and 'help timap extended'.
  * \author Nathan D. Levinzon <ndlevinzon@gmail.com>
  */
class Exec_TIMap : public Exec {
  public:
    Exec_TIMap() : Exec(PARM) {}
    void Help() const;
    void Help(ArgList&) const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_TIMap(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    /// Aligned-structure output format (default Amber OFF .lib for LEaP/TI).
    enum AlignedOutFmt { AOUT_LIB = 0, AOUT_MOL2, AOUT_PDB };
    /// One topology in a series-mode parent search.
    struct SeriesMember {
      std::string name;
      Topology* top;
      DataSet* ds;
    };

    Topology* FindNamedTop(DataSetList&, std::string const&, DataSet**) const;
    void ParseTopSpec(std::string const&, std::string&, std::string&) const;
    std::string FileStem(FileName const&) const;
    Topology* ResolveTop(CpptrajState&, std::string const&, DataSet**) const;
    bool LoadCoords(Topology const&, DataSet*, Frame&) const;
    void CopyXyz(Frame const&, int, Frame&, int) const;
    void ApplyOrderToFrame(Frame const&, std::vector<int> const&, Frame&) const;
    std::string OffUnitName(std::string, char) const;
    NameType ResNameOf(Topology const&, const char*) const;
    int BuildTiUnit(bool, Topology const&, Topology const&, Frame const&, Frame const&,
                    std::vector<TIMatch::Result::DualSlot> const&,
                    Topology&, Frame&, NameType const&) const;
    int WriteMol2File(std::string const&, Topology&, Frame const&, DataSetList const&) const;
    int WritePdbFile(std::string const&, Topology&, Frame const&, DataSetList const&) const;
    const char* AlignedOutExt(AlignedOutFmt) const;
    const char* AlignedOutLabel(AlignedOutFmt) const;
    int ParseAlignedOutFmt(std::string const&, AlignedOutFmt&) const;
    AlignedOutFmt FmtFromFilename(std::string const&, AlignedOutFmt) const;
    std::string Q(std::string const&) const;
    std::string AtomTypeStr(Atom const&) const;
    int AtomNum1(Topology const&, const char*) const;
    void LeapConnect(Topology const&, int&, int&, const char*&) const;
    int WriteAmberLib(std::string const&, std::string const&, Topology const&, Frame const&) const;
    int WriteAlignedOut(std::string const&, AlignedOutFmt, Topology&, Frame const&,
                        DataSetList const&) const;
    int WriteScmask(std::string const&, Topology const&, Topology const&,
                    std::vector<TIMatch::Result::DualSlot> const&) const;
    int WriteDualAtoms(std::string const&, Topology const&, Topology const&,
                       std::vector<TIMatch::Result::DualSlot> const&) const;
    int WriteMapFile(std::string const&, Topology const&, Topology const&,
                     TIMatch::Result const&, std::string const&, std::string const&) const;
    std::string SafeFileToken(std::string const&) const;
    void PushSeriesName(std::vector<std::string>&, std::string const&) const;
    double LeadoptSimScore(int, int) const;
    int SelectSeriesParent(TIMatch const&, std::vector<SeriesMember> const&,
                           std::vector<double>&, std::vector<int>&) const;
};
#endif
