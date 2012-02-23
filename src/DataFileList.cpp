// DataFileList
#include <cstddef> // NULL
#include "DataFileList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataFileList::DataFileList() {
//  fprintf(stderr,"DataFileList CONSTRUCTOR\n");
  debug=0;
}

// DESTRUCTOR
DataFileList::~DataFileList() {
//  fprintf(stderr,"DataFileList DESTRUCTOR\n");
  for (it = fileList.begin(); it != fileList.end(); it++)
    delete *it;
}

// DataFileList::SetDebug()
/** Set debug level for DataFileList and all datafiles in it. */
void DataFileList::SetDebug(int debugIn) {
  debug=debugIn;
  if (debug>0)
    mprintf("DataFileList DEBUG LEVEL SET TO %i\n",debug);
}

// DataFileList::GetDataFile()
/** Return DataFile specified by given file name if it exists in the list,
  * otherwise return null.
  */
DataFile *DataFileList::GetDataFile(char *nameIn) {
  DataFile *Current;

  if (nameIn==NULL) return NULL;
  Current=NULL;
  for (it = fileList.begin(); it != fileList.end(); it++) {
    if ( (*it)->DataFileNameIs(nameIn) ) {
      Current = *it;
      break;
    }
  }
  
  return Current;
}

// DataFileList::Add()
/** Add dataset to datafile in list with given file name. If the file does
  * not yet exist in the list create it. Return a pointer to the datafile
  * in the list.
  */
DataFile *DataFileList::Add(char *nameIn, DataSet *D) {
  DataFile *Current;
  //char tempName[1024]; // DEBUG

  // If no filename, no output desired
  if (nameIn==NULL) return NULL;
  // If DataSet is NULL, dont add
  if (D==NULL) {
    mprintf("Error: Attempting to add non-existent dataset to file %s\n",nameIn);
    return NULL;
  }

  // Append thread prefix to filename
  //sprintf(tempName,"%s.%03i",nameIn,worldrank); // DEBUG

  // Check if this filename already in use
  Current = this->GetDataFile(nameIn);

  // If no DataFile associated with nameIn, create new datafile
  if (Current==NULL) {
    Current = new DataFile(nameIn); 
    fileList.push_back(Current);
  }

  // Add the dataset to the current DataFile
  Current->AddSet(D);

  // Set debug level
  Current->SetDebug(debug);

  // DEBUG
  //mprintf("** ADDED DATASET %s TO FILE %s\n",D->Name(),Current->filename);

  return Current;
}

// DataFileList::Info()
/** Print information on what datasets are going to what datafiles */
void DataFileList::Info() {

  if (fileList.empty()) {
    //mprintf("NO DATASETS WILL BE OUTPUT\n");
    return;
  }

  mprintf("DATAFILE OUTPUT:\n");
  for (it = fileList.begin(); it != fileList.end(); it++) {
    mprintf("  %s: ",(*it)->Filename());
    (*it)->DataSetNames();
    mprintf("\n");
  }
}

// DataFileList::Write()
/** Call write for all datafiles in list. Only master should call this.
  */
void DataFileList::Write() {

  //if (worldrank!=0) return; 

  for (it = fileList.begin(); it != fileList.end(); it++)
    (*it)->Write();
}

// DataFileList::ProcessDataFileArgs()
/** Process command relating to data files. All datafile commands have format:
  *   datafile <cmd> <datafile> ...
  */
void DataFileList::ProcessDataFileArgs(DataSetList *masterDSL) {
  char *df_cmd = NULL;
  char *name1 = NULL;
  char *name2 = NULL;
  int width,precision;
  DataFile *df;

  if (DF_Args.empty()) return;
  mprintf("DATAFILE SETUP:\n");

  // Loop through all "datafile" arguments
  for (std::list<ArgList>::iterator dataArg=DF_Args.begin();
                                    dataArg!=DF_Args.end();
                                    dataArg++)
  {
    // Next string will be the argument passed to datafile
    df_cmd = (*dataArg).getNextString();
    mprintf("  [%s]\n",(*dataArg).ArgLine());
    // Next string is datafile that command pertains to.
    name1 = (*dataArg).getNextString();
    if (name1==NULL) {
      mprintf("Error: datafile %s: No filename given.\n",df_cmd);
      continue;
    }
    df = GetDataFile(name1);

    // datafile create
    // Usage: datafile create <filename> <dataset0> <dataset1> ...
    if ( (*dataArg).ArgIs(1,"create") ) {
      if (df==NULL)
        mprintf("    Creating file %s\n",name1);
      while ( (name2=(*dataArg).getNextString())!=NULL ) {
        if ( Add(name1, masterDSL->Get(name2))==NULL ) {
          mprintf("Warning: Dataset %s does not exist in main dataset list, skipping.\n",name2);
        }
      }

    // datafile xlabel
    // Usage: datafile xlabel <filename> <label>
    } else if ( (*dataArg).ArgIs(1,"xlabel") ) {
      if (df==NULL) {
        mprintf("Error: datafile xlabel: DataFile %s does not exist.\n",name1);
        continue;
      }
      df->SetXlabel((*dataArg).getNextString());

    // datafile ylabel
    // Usage: datafile ylabel <filename> <label>
    } else if ( (*dataArg).ArgIs(1,"ylabel") ) {
      if (df==NULL) {
        mprintf("Error: datafile ylabel: DataFile %s does not exist.\n",name1);
        continue;
      }
      df->SetYlabel((*dataArg).getNextString());

    // datafile time
    // Usage: datafile time <datafile> <time>
    } else if ( (*dataArg).ArgIs(1,"time") ) {
      if (df==NULL) {
        mprintf("Error: datafile time: DataFile %s does not exist.\n",name1);
        continue;
      }
      df->SetXstep((*dataArg).getNextDouble(1));

    // datafile invert
    // Usage: datafile invert <filename>
    } else if ( (*dataArg).ArgIs(1,"invert") ) {
      if (df==NULL) {
        mprintf("Error: datafile invert: DataFile %s does not exist.\n",name1);
        continue;
      }
      mprintf("    Inverting datafile %s\n",name1);
      df->SetInverted();

    // datafile noxcol
    // Usage: datafile noxcol <filename>
    } else if ( (*dataArg).ArgIs(1,"noxcol") ) {
      if (df==NULL) {
        mprintf("Error: datafile noxcol: DataFile %s does not exist.\n",name1);
        continue;
      }
      mprintf("    Not printing x column for datafile %s\n",name1);
      df->SetNoXcol();

    // datafile noheader
    // Usage: datafile noheader <filename>
    } else if ( (*dataArg).ArgIs(1,"noheader") ) {
      if (df==NULL) {
        mprintf("Error: datafile noheader: DataFile %s does not exist.\n",name1);
        continue;
      }
      mprintf("    Not printing header for datafile %s\n",name1);
      df->SetNoHeader();

    // datafile precision
    // Usage: datafile precision <filename> <dataset> [<width>] [<precision>]
    //        If width/precision not specified default to 12.4
    } else if ( (*dataArg).ArgIs(1,"precision") ) {
      if (df==NULL) {
        mprintf("Error: datafile precision: DataFile %s does not exist.\n",name1);
        continue;
      }
      // This will break if dataset name starts with a digit...
      width = (*dataArg).getNextInteger(12);
      precision = (*dataArg).getNextInteger(4);
      name2 = (*dataArg).getNextString();
      df->SetPrecision(name2,width,precision);
    }

  } // END loop over datafile args
}

