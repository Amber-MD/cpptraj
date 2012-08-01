#include "DataIO_Grace.h"

// CONSTRUCTOR
DataIO_Grace::DataIO_Grace() {
  ymin_ = 1;
  ystep_ = 1;
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

