#ifndef INC_DATAIO_AMBERLIB_H
#define INC_DATAIO_AMBERLIB_H
#include "DataIO.h"
class BufferedLine;
/// <Enter description of DataIO_AmberLib here>
class DataIO_AmberLib : public DataIO {
  public:
    DataIO_AmberLib();
    static void ReadHelp();
    static void WriteHelp();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_AmberLib(); }
    int processReadArgs(ArgList&);
    int ReadData(FileName const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(FileName const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&);
  private:
    // Keep in sync with sectionStr_
    enum SectionType { ATOMTABLE = 0, PERTINFO, BOUNDBOX, CHILDSEQUENCE, CONNECT,
                       CONNECTIVITY, HIERARCHY, UNITNAME, POSITIONS, RESCONNECT,
                       RESIDUES, SOLVENTCAP, VELOCITIES, UNKNOWN_SECTION };
    /// Strings corresponding to section type
    static const char* sectionStr_[];

    /// ID OFF section from line
    static inline SectionType id_section(std::string const&, std::string const&);
    /// Read unit from OFF file
    int read_unit(BufferedLine&, std::string&, std::string const&) const;
};
#endif
