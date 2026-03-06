#ifndef INC_DATASET_PDBRESMAP_H
#define INC_DATASET_PDBRESMAP_H
#include "DataSet.h"
#include <map>
#include <vector>
#include <string>
#include "PdbResMapType.h"
/// Hold LEaP PDB residue name map 
class DataSet_PdbResMap : public DataSet {
    typedef std::vector<std::string> Sarray;
    typedef std::pair<NameType, Sarray> PairType;
    typedef std::map<NameType, Sarray> MapType;
  public:
    DataSet_PdbResMap();
    static DataSet* Alloc() { return (DataSet*)new DataSet_PdbResMap(); }
    // ----- DataSet functions -------------------
    size_t Size()                                    const { return 0; }
    void Info()                                      const { return; }
    int Allocate(SizeArray const&)                         { return 1; }
    void Add(size_t, const void*)                          { return; }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
    int Append(DataSet*)                                   { return 1; }
    size_t MemUsageInBytes()                         const { return 0; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    // -------------------------------------------
    /// Add a PDB name to unit name mapping
    int AddPdbResMap(Cpptraj::PdbResMapType const&);
    /// Find unit name based on PDB name/terminal type
    std::string FindUnitName(NameType const&, Cpptraj::Structure::TerminalType) const;

    void PrintPdbResMap() const;
  private:
    MapType pdbResMap_; ///< Map PDB residue name to array of unit names: begin, non-terminal, end
};
#endif
