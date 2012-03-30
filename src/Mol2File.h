#ifndef INC_MOL2FILE_H
#define INC_MOL2FILE_H
#include "Topology.h"
#include "FileIO.h"
/// Used to access mol2 files.
class Mol2File {
  public: 
    Mol2File();

    enum TRIPOSTAG { MOLECULE=0, ATOM, BOND, SUBSTRUCT };

    bool IsMol2Keyword();
    //bool NextLine(FileIO*);
    int ScanTo( FileIO *, TRIPOSTAG );
    Atom Mol2Atom();
    Residue Mol2Residue();
  private:
    static const char TRIPOSTAGTEXT[][22];
  protected: // TODO: Should be private
    static const size_t BUF_SIZE_ = 256;
    char buffer_[BUF_SIZE_];
};
#endif  
