#ifndef INC_DATAFILELIST_H
#define INC_DATAFILELIST_H
#include "DataFile.h"
#include "DataSet.h"
#include "ArgList.h"
/// Holds a list of output DataFiles/CpptrajFiles.
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
    void SetEnsembleNum(int i) { ensembleNum_ = i; }
    void SetEnsExtension(bool b) { ensExt_ = b; }
    /// \return DataFile whose full path matches given string or 0.
    DataFile* GetDataFile(FileName const&) const;
    /// \return CpptrajFile whose full path matches given string or 0.
    CpptrajFile* GetCpptrajFile(FileName const&) const;
    /// \return DataFile specified by name, add if none exists, or 0 if no name specified.
    DataFile* AddDataFile(FileName const&, ArgList&);
    /// \return DataFile specified by name, add if none exists, or 0 if no name specified.
    DataFile* AddDataFile(FileName const&);
    /// \return DataFile specified by name with specific format, add if none exists.
    DataFile* AddDataFile(FileName const&, ArgList&, DataFile::DataFormatType);
    /// \return DataFile specified by name with specific format, add if none exists.
    DataFile* AddDataFile(FileName const&, DataFile::DataFormatType, ArgList const&);
    /// Types of CpptrajFile that can be created.
    enum CFtype { TEXT = 0, PDB };
    /// \return CpptrajFile specified by name, add if none exists, or 0 if no name specified.
    CpptrajFile* AddCpptrajFile(FileName const&,std::string const&);
    /// Add/create CpptrajFile of given type. No STDOUT.
    CpptrajFile* AddCpptrajFile(FileName const&,std::string const&,CFtype);
    /// Add/create CpptrajFile of given type; optionally allow STDOUT.
    CpptrajFile* AddCpptrajFile(FileName const&,std::string const&,CFtype,bool);
    /// List DataFiles and CpptrajFiles.
    void List() const;
    /// Write all DataFiles in list that have not yet been written.
    void WriteAllDF();
    /// \return true if DataFiles have not yet been written.
    bool UnwrittenData() const;
    void ResetWriteStatus();
    int ProcessDataFileArgs(ArgList&);
    int Debug() const { return debug_; }
    int EnsembleNum() const { return ensembleNum_; }
  private:
    int GetCpptrajFileIdx(FileName const&) const;

    typedef std::vector<DataFile*> DFarray;
    DFarray fileList_;
    // NOTE: CpptrajFile* must be kept in its own array since currently the
    //       copy constructor copies them closed no matter what.
    typedef std::vector<CpptrajFile*> CFarray;
    CFarray cfList_;
    /// Class to hold CpptrajFile metadata.
    class CFstruct;
    typedef std::vector<CFstruct> MDarray;
    MDarray cfData_;

    int debug_;
    int ensembleNum_; ///< Ensemble member number.
    bool ensExt_;     ///< If true append ensemble member number to file names.
};
// ----- INTERNAL CLASS DEFINITIONS --------------------------------------------
class DataFileList::CFstruct {
  public:
    CFstruct() : type_(TEXT) {}
    CFstruct(std::string const& d, CFtype t) : description_(d), type_(t) {}
    const char* descrip()        const { return description_.c_str(); }
    CFtype Type()                const { return type_;                }
    void UpdateDescrip(std::string const& d) { description_.append(", " + d); }
  private:
    std::string description_; ///< Description of file contents.
    CFtype type_;             ///< The file subtype.
};
#endif
