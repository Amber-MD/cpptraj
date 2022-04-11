#ifndef INC_FORLOOP_DATASETBLOCKS_H
#define INC_FORLOOP_DATASETBLOCKS_H
#include "ForLoop.h"
/// For loop, data set blocks
class ForLoop_dataSetBlocks : public ForLoop {
  public:
    ForLoop_dataSetBlocks();

    static void helpText();

    int SetupFor(CpptrajState&, ArgList&);
    int BeginFor(DataSetList const&);
    // NOTE: Not a const ref so dataSetBlocks loop can create sets
    bool EndFor(DataSetList&);
    
  private:
    enum ModeType { BLOCKS=0, CUMULATIVE };

    DataSet* sourceSet_;   ///< Set to loop over
    DataSet* currentSet_;  ///< Current subset
    std::string sourceSetName_; 
    long int blocksize_;   ///< Size of blocks
    long int blockoffset_; ///< Block offset
    long int idx_;         ///< Current index into sourceSet_
    ModeType mode_;        ///< Current mode
};
#endif
