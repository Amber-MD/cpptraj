#include "DataIO_Grace.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataIO_Grace::DataIO_Grace() {
  ymin_ = 1;
  ystep_ = 1;
}

// Dont assume anything about set ordering
int DataIO_Grace::ReadData(DataSetList& datasetlist) {
  ArgList dataline;
  char buffer[1024];
  int setnum = 0;
  DataSet *dset = NULL;
  std::vector<DataSet*> Dsets;
  std::vector<std::string> labels;

  while ( IO->Gets(buffer, 1024) == 0 ) {
    if ( buffer[0] == '@' ) { // FIXME: Dont assume it is the first char
      // Create command line without the @
      dataline.SetList(buffer+1, " \t");
      if ( dataline[1] == "legend" ) {
        labels.push_back( dataline[2] );
      } else if (dataline[0] == "target") {
        dset = datasetlist.AddSetIdx( DataSet::DOUBLE, BaseFileName(), setnum++);
        if (dset == NULL) {
          mprinterr("Error: %s: Could not allocate data set.\n", Name());
          return 1;
        }
        Dsets.push_back( dset );
      }
    } else {
      // Data
      if (dset==NULL) {
        mprinterr("Error: %s: Malformed grace file. Expected 'target' before data.\n", Name());
        return 1;
      }
      dataline.SetList(buffer, " \t");
      int frame = convertToInteger( dataline[0] );
      double dval = convertToDouble( dataline[1] );
      dset->Add( frame, &dval );
    }
  }

  if (!labels.empty()) {
    if (Dsets.size() == labels.size()) {
      mprintf("\tLabels:\n");
      std::vector<DataSet*>::iterator set = Dsets.begin();
      for (std::vector<std::string>::iterator lbl = labels.begin();
                                              lbl != labels.end(); ++lbl)
      {
        mprintf("\t\t%s\n", (*lbl).c_str());
        (*set)->SetLegend( *lbl );
        ++set;
      }
    } else {
      mprintf("Warning: Number of labels (%zu) does not match number of sets (%zu)\n",
              labels.size(), Dsets.size());
    }
  }
    
  return 0;
}
      
// DataIO_Grace::processWriteArgs()
int DataIO_Grace::processWriteArgs(ArgList &argIn) {
  char *ylabel = argIn.getKeyString("ylabel",NULL);
  if (ylabel!=NULL) y_label_.assign(ylabel);
  ymin_ = argIn.getKeyDouble("ymin",ymin_);
  ystep_ = argIn.getKeyDouble("ystep",ystep_);
  return 0;
}

// DataIO_Grace::WriteData()
int DataIO_Grace::WriteData(DataSetList &SetList) {
  DataSetList::const_iterator set;
  // Create format string for X column. Default precision is 3
  SetupXcolumn();

  // Grace header
  Printf(
    "@with g0\n@  xaxis label \"%s\"\n@  yaxis label \"%s\"\n@  legend 0.2, 0.995\n@  legend char size 0.60\n",
    x_label_.c_str(), y_label_.c_str()
  );
  // Loop over sets in data
  int setnum = 0;
  for (set=SetList.begin(); set != SetList.end(); set++) {
    // Set information
    Printf("@  s%-8i legend \"%s\"\n@target G0.S%-8i\n@type xy\n", 
                   setnum, (*set)->Legend().c_str(), setnum );
    // Write Data
    for (int frame=0; frame<maxFrames_; frame++) {
      // If specified, run through every set in the frame and check if empty
      if (!printEmptyFrames_) {
        if ( (*set)->FrameIsEmpty(frame) ) continue;
      }
      double xcoord = (xstep_ * (double)frame) + xmin_;
      Printf(x_format_.c_str(), xcoord);
      (*set)->WriteBuffer(*this,frame);
      Printf("\n");
    }
    ++setnum;
  }
  return 0;
}

// DataIO_Grace::WriteDataInverted()
int DataIO_Grace::WriteDataInverted(DataSetList &SetList) {
  DataSetList::const_iterator set;
  // Create format string for X column. Default precision is 3
  SetupXcolumn();

  // Grace Header
  Printf(
    "@with g0\n@  xaxis label \"%s\"\n@  yaxis label \"%s\"\n@  legend 0.2, 0.995\n@  legend char size 0.60\n",
    y_label_.c_str(), x_label_.c_str()
  );

  // Loop over frames
  for (int frame=0; frame<maxFrames_; frame++) {
    // If specified, run through every set in the frame and check if empty
    if (!printEmptyFrames_) {
      bool hasEmpty=false;
      for (set=SetList.begin(); set!=SetList.end(); set++) {
        if ( (*set)->FrameIsEmpty(frame) ) {
          hasEmpty=true; 
          break;
        }
      }
      if (hasEmpty) continue;
    }
    // Set information
    Printf("@target G0.S%-8i\n@type xy\n",frame);
    // Loop over all Set Data for this frame
    int setnum = 0;
    for (set=SetList.begin(); set!=SetList.end(); set++) {
      double xcoord = (ystep_ * setnum) + ymin_;
      Printf(x_format_.c_str(), xcoord);
      (*set)->WriteBuffer(*this, frame);
      Printf("\"%s\"\n", (*set)->c_str());
    }
  }
  return 0;
}

