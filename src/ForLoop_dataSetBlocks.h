#ifndef INC_FORLOOP_DATASETBLOCKS_H
#define INC_FORLOOP_DATASETBLOCKS_H
#include "ForLoop.h"
/// For loop, data set blocks
class ForLoop_dataSetBlocks : public ForLoop {
  public:
    ForLoop_dataSetBlocks();

    int SetupFor(CpptrajState&, ArgList&);
    int BeginFor(DataSetList const&);
    bool EndFor(DataSetList const&);
    
  private:
    DataSet* sourceSet_;   ///< Set to loop over
    DataSet* currentSet_;  ///< Current subset
    std::string sourceSetName_; 
    long int blocksize_;   ///< Size of blocks
    long int blockoffset_; ///< Block offset
    long int idx_;         ///< Current index into sourceSet_
};
#endif
