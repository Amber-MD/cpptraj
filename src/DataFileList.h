#ifndef INC_DATAFILELIST_H
#define INC_DATAFILELIST_H
#include "DispatchObject.h"
#include "FileList.h"
#include "DataFile.h"
#include "DataSet.h"
#include "ArgList.h"
// Class: DataFileList
/// Holds a list of DataFile classes. 
class DataFileList : public FileList {
  public:
    DataFileList();
    ~DataFileList();
    void Clear();
    void SetDebug(int);
#   ifdef MPI
    void SetEnsembleMode(int mIn) { ensembleMode_ = mIn; }
#   endif
    DataFile* GetDataFile(std::string const&) const;
    DataFile* AddDataFile(std::string const&, ArgList&);
    DataFile* AddDataFile(std::string const&);
    // TODO: Deprecate in favor of AddDataFile
    DataFile* AddSetToFile(std::string const&,  DataSet*);
    void List() const;
    /// Write all DataFiles in list that have not yet been written.
    void WriteAllDF();
    int ProcessDataFileArgs(ArgList&);
  private:
    typedef std::vector<DataFile*> DFarray;
    DFarray fileList_;
    int debug_;
#   ifdef MPI
    int ensembleMode_; ///< When parallel reading, append filenames with this
#   endif
};
#endif
