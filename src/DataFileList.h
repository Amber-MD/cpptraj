#ifndef INC_DATAFILELIST_H
#define INC_DATAFILELIST_H
#include "DataFile.h"
#include "DataSet.h"
#include "ArgList.h"
// Class: DataFileList
/// Holds a list of DataFile classes.
/** The DataFileList is meant to hold all output data files defined by any
  * Actions or Analysis. This allows multiple sets to be directed to the
  * same file etc. Holds both DataFiles, which hold DataSets, and CpptrajFiles
  * which are intended for plain ASCII text output (like from 'hbond avgout').
  */ 
class DataFileList {
  public:
    DataFileList();
    ~DataFileList();
    void Clear();
    DataFile* RemoveDataFile(DataFile*);
    void RemoveDataSet(DataSet*);
    void SetDebug(int);
#   ifdef MPI
    void SetEnsembleMode(int);
#   endif
    /// \return DataFile whose full path matches given string or 0.
    DataFile* GetDataFile(std::string const&) const;
    /// \return CpptrajFile whose full path matches given string or 0.
    CpptrajFile* GetCpptrajFile(std::string const&) const;
    /// Add DataFile to list if name specified, or return already existing DataFile.
    DataFile* AddDataFile(std::string const&, ArgList&);
    /// Add DataFile to list if name specified, or return already existing DataFile.
    DataFile* AddDataFile(std::string const&);
    // TODO: Deprecate in favor of AddDataFile
    DataFile* AddSetToFile(std::string const&,  DataSet*);
    /// Types of cpptrajfile that can be created.
    enum CFtype { TEXT = 0, PDB };
    /// Add CpptrajFile to list if name specified, or return already existing CpptrajFile.
    CpptrajFile* AddCpptrajFile(std::string const&,std::string const&);
    CpptrajFile* AddCpptrajFile(std::string const&,std::string const&,CFtype);
    CpptrajFile* AddCpptrajFile(std::string const&,std::string const&,CFtype,bool);
    /// Open any pending CpptrajFiles.
    int OpenCpptrajFiles();
    /// Close any open CpptrajFiles.
    void CloseCpptrajFiles();
    /// List DataFiles and CpptrajFiles.
    void List() const;
    /// Write all DataFiles in list that have not yet been written.
    void WriteAllDF();
    void ResetWriteStatus();
    int ProcessDataFileArgs(ArgList&);
    int Debug() const { return debug_; }
  private:
    int GetCpptrajFileIdx(std::string const&) const;

    typedef std::vector<DataFile*> DFarray;
    DFarray fileList_;
    // NOTE: CpptrajFile* must be kept in its own array since currently the
    //       copy constructor copies them closed no matter what.
    typedef std::vector<CpptrajFile*> CFarray;
    CFarray cfList_;
    /** FIRSTOPEN: File should be opened write.
      * REOPEN: File should be opened append.
      * NO_ACTION: File does not need to be opened.
      */
    enum CFmode { FIRSTOPEN = 0, REOPEN, NO_ACTION };
    ///< Class to hold CpptrajFile metadata.
    class CFstruct;
    typedef std::vector<CFstruct> MDarray;
    MDarray cfData_;

    int debug_;
};
// ----- INTERNAL CLASS DEFINITIONS --------------------------------------------
class DataFileList::CFstruct {
  public:
    CFstruct() : type_(TEXT), mode_(NO_ACTION) {}
    CFstruct(std::string const& d, CFtype t) : description_(d), type_(t), mode_(FIRSTOPEN) {}
    const char* descrip()        const { return description_.c_str(); }
    CFtype Type()                const { return type_;                }
    CFmode Mode()                const { return mode_;                }
    void UpdateDescrip(std::string const& d) { description_.append(", " + d); }
    void UpdateStatus( CFmode m )            { mode_ = m;                     }
  private:
    std::string description_; ///< Description of file contents.
    CFtype type_;             ///< The file subtype.
    CFmode mode_;             ///< What to do with file on next open.
};
#endif
