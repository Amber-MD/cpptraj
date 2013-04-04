#include <cstdio> // sscanf
#include "DataIO_Grace.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"

// CONSTRUCTOR
DataIO_Grace::DataIO_Grace() :
  ymin_(1),
  ystep_(1)
{}

// Dont assume anything about set ordering
int DataIO_Grace::ReadData(std::string const& fname, DataSetList& datasetlist) {
  ArgList dataline;
  int setnum = 0;
  int frame = 0;
  DataSet *dset = 0;
  std::vector<DataSet*> Dsets;
  std::vector<std::string> labels;
  double dval;
  const char* linebuffer;
  
  // Allocate and set up read buffer
  BufferedLine buffer;
  if (buffer.OpenRead( fname )) return 1;
  buffer.SetupBuffer();

  // Read chunks from file
  while ( (linebuffer = buffer.Line()) != 0 ) {
    if (linebuffer[0] == '@') {
      // Command: create command line without the @
      dataline.SetList(linebuffer+1, " \t");
      if ( !dataline.CommandIs("legend") && dataline.Contains("legend") ) {
        // Legend keyword that does not come first. Dont store blanks.
        std::string lbl = dataline.GetStringKey("legend");
        if (!lbl.empty())
          labels.push_back( lbl );
      } else if (dataline.CommandIs("target")) {
        // Indicates dataset will be read soon. Allocate new set.
        dset = datasetlist.AddSetIdx( DataSet::DOUBLE, buffer.Filename().Base(), setnum++);
        if (dset == 0) {
          mprinterr("Error: %s: Could not allocate data set.\n", buffer.Filename().full());
          return 1;
        }
        Dsets.push_back( dset );
        frame = 0;
      }
    } else if (linebuffer[0] != '#' && linebuffer[0] != '&') { // Skip comments and dataset end
      // Data
      if (dset==0) {
        mprinterr("Error: %s: Malformed grace file. Expected 'target' before data.\n", 
                  buffer.Filename().full());
        return 1;
      }
      // FIXME: Ignoring frame for now
      sscanf(linebuffer,"%*s %lf", &dval);
      dset->Add( frame++, &dval );
    }
  } // END loop over file
  buffer.CloseFile();

  // Set DataSet legends if specified
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
  std::string ylabel = argIn.GetStringKey("ylabel");
  if (!ylabel.empty()) y_label_.assign(ylabel);
  ymin_ = argIn.getKeyDouble("ymin",ymin_);
  ystep_ = argIn.getKeyDouble("ystep",ystep_);
  return 0;
}

// DataIO_Grace::WriteData()
int DataIO_Grace::WriteData(std::string const& fname, DataSetList &SetList) {
  DataSetList::const_iterator set;
  CpptrajFile file;
  if (file.OpenWrite(fname)) return 1;
  // Create format string for X column. Default precision is 3
  SetupXcolumn();

  // Grace header
  file.Printf(
    "@with g0\n@  xaxis label \"%s\"\n@  yaxis label \"%s\"\n@  legend 0.2, 0.995\n@  legend char size 0.60\n",
    x_label_.c_str(), y_label_.c_str()
  );
  // Loop over sets in data
  int setnum = 0;
  for (set=SetList.begin(); set != SetList.end(); set++) {
    // Set information
    file.Printf("@  s%-8i legend \"%s\"\n@target G0.S%-8i\n@type xy\n", 
                   setnum, (*set)->Legend().c_str(), setnum );
    // Write Data
    for (int frame=0; frame<maxFrames_; frame++) {
      // If specified, run through every set in the frame and check if empty
      if (!printEmptyFrames_) {
        if ( (*set)->FrameIsEmpty(frame) ) continue;
      }
      PrintX( file, frame );
      (*set)->WriteBuffer(file, frame);
      file.Printf("\n");
    }
    ++setnum;
  }
  file.CloseFile();
  return 0;
}

// DataIO_Grace::WriteDataInverted()
int DataIO_Grace::WriteDataInverted(std::string const& fname, DataSetList &SetList) {
  DataSetList::const_iterator set;
  CpptrajFile file;
  if (file.OpenWrite(fname)) return 1;
  // Create format string for X column. Default precision is 3
  SetupXcolumn();

  // Grace Header
  file.Printf(
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
    file.Printf("@target G0.S%-8i\n@type xy\n",frame);
    // Loop over all Set Data for this frame
    int setnum = 0;
    for (set=SetList.begin(); set!=SetList.end(); set++) {
      double xcoord = (ystep_ * setnum) + ymin_;
      file.Printf(x_format_.c_str(), xcoord);
      (*set)->WriteBuffer(file, frame);
      file.Printf("\"%s\"\n", (*set)->Legend().c_str());
    }
  }
  file.CloseFile();
  return 0;
}

