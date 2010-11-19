#ifndef INC_DATASETLIST_H
#define INC_DATASETLIST_H
// DataSetList
/* Data sets are added to the list by various actions. A data file is created
 * and the set is added to the data file list if an out keyword was specified
 * in the action argument line.
 */

#include "DataSet.h"
//#include "DataFileList.h"

class DataSetList {

    DataSet **DataList;
    int Ndata;

  public:
    int maxFrames; // The maximum number of frames to be read, set in PtrajState::Run

    DataSetList();
    ~DataSetList();

    DataSet *Get(char *);
//    DataSet *Add( dataType , ArgList*);
    DataSet *Add( dataType , char*);
    int AddData(int , void *, int );
    void Info();
    void Sync();
};

#endif
