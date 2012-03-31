#include "DataIO_Gnuplot.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
GnuplotDataFile::GnuplotDataFile() {
  y_label_ = "";
  ymin_ = 1;
  ystep_ = 1;
  useMap_ = false;
  printLabels_ = true;
}

// GnuplotDataFile::processWriteArgs()
int GnuplotDataFile::processWriteArgs(ArgList &argIn) {
  char *ylabel = argIn.getKeyString("ylabel",NULL);
  if (ylabel!=NULL) y_label_.assign(ylabel);
  ymin_ = argIn.getKeyDouble("ymin",ymin_);
  ystep_ = argIn.getKeyDouble("ystep",ystep_);

  if (argIn.hasKey("nolabels")) printLabels_ = false;
  if (argIn.hasKey("usemap")) useMap_ = true;
  return 0;
}

// GnuplotDataFile::WriteData()
/** Write each frame from all sets in blocks in the following format:
  *   Frame Set   Value
  * Originally there was a -0.5 offset for the Set values in order to center
  * grid lines on the Y values, e.g.
  *   1     0.5   X
  *   1     1.5   X
  *   1     2.5   X
  *
  *   2     0.5   X
  *   2     1.5   X
  *   ...
  * However, in the interest of keeping data consistent, this is no longer
  * done. Could be added back in later as an option.
  */
int GnuplotDataFile::WriteData(DataSetList &SetList) {
  CharBuffer buffer;
  DataSetList::const_iterator set;
  double xcoord, ycoord;

  // Create format string for X and Y columns. Default precision is 3
  SetupXcolumn();
  std::string xy_format_string = x_format_ + " " + x_format_ + " ";
  const char *xy_format = xy_format_string.c_str();

  // Turn off labels if number of sets is too large since they 
  // become unreadable. Should eventually have some sort of 
  // autotick option.
  if (printLabels_ && SetList.Size() > 30 ) {
    mprintf("Warning: %s: gnuplot - number of sets %i > 30, turning off Y labels.\n",
            BaseName(), SetList.Size());
    printLabels_ = false;
  }

  // Map command
  if (!useMap_)
    Printf("set pm3d map corners2color c1\n");
  else
    Printf("set pm3d map\n");

  // Y axis Data Labels
  if (printLabels_) {
    // NOTE: Add option to turn on grid later?
    //outfile->Printf("set pm3d map hidden3d 100 corners2color c1\n");
    //outfile->Printf("set style line 100 lt 2 lw 0.5\n");
    // Set up Y labels
    Printf("set ytics %8.3f,%8.3f\nset ytics(",ymin_,ystep_);
    int setnum = 0;
    for (set=SetList.begin(); set!=SetList.end(); set++) {
      if (setnum>0) Printf(",");
      ycoord = (ystep_ * setnum) + ymin_;
      Printf("\"%s\" %8.3f",(*set)->Name(),ycoord);
      ++setnum;
    }
    Printf(")\n");
  }

  // Set axis label and range, write plot command
  // Make Yrange +1 and -1 so entire grid can be seen
  ycoord = (ystep_ * SetList.Size()) + ymin_;
  // Make Xrange +1 and -1 as well
  xcoord = (xstep_ * maxFrames_) + xmin_;
  Printf
    ("set xlabel \"%s\"\nset ylabel \"%s\"\nset yrange [%8.3f:%8.3f]\nset xrange [%8.3f:%8.3f]\nsplot \"-\" with pm3d title \"%s\"\n",
    x_label_.c_str(), y_label_.c_str(),
    ymin_ - ystep_, ycoord + ystep_,
    xmin_ - xstep_, xcoord + xstep_, BaseName()
  );

  // Allocate space for data.
  // Each frame: (Nset * ((XY + setwidth + 1) [+ (XY + 2)]:!useMap) + 1 (newline)
  size_t dataFileSize = 0;
  for (set=SetList.begin(); set !=SetList.end(); set++) {
    size_t dataSetSize = (xcol_width_ + xcol_width_ + 3 + (*set)->Width());
    if (!useMap_)
      dataSetSize += (xcol_width_ + xcol_width_ + 4);
    ++dataSetSize; // Newline
    dataFileSize += ((dataSetSize * maxFrames_) + 1);
  }
  // If !useMap, += ((Nset * (XY + 2) + 1)
  if (!useMap_)
    dataFileSize += ((SetList.Size() * (xcol_width_ + xcol_width_ + 4))+1);
  // + 1 (terminal newline)
  ++dataFileSize;
  buffer.Allocate( dataFileSize );

  // Data
  int frame = 0;
  for (; frame < maxFrames_; frame++) {
    xcoord = (xstep_ * frame) + xmin_;
    int setnum = 0;
    for (set=SetList.begin(); set !=SetList.end(); set++) {
      ycoord = (ystep_ * setnum) + ymin_;
      buffer.WriteXY(xy_format, xcoord, ycoord);
      (*set)->WriteBuffer(buffer,frame);
      buffer.NewLine();
      ++setnum;
    }
    if (!useMap_) {
      // Print one empty row for gnuplot pm3d without map
      ycoord = (ystep_ * setnum) + ymin_;
      buffer.WriteXY(xy_format, xcoord, ycoord);
      buffer.WriteInteger("%1i",0);
      buffer.NewLine();
    }
    buffer.NewLine();
  }
  if (!useMap_) {
    // Print one empty set for gnuplot pm3d without map
    xcoord = (xstep_ * frame) + xmin_;
    for (int blankset=0; blankset <= SetList.Size(); blankset++) {
      ycoord = (ystep_ * blankset) + ymin_;
      buffer.WriteXY(xy_format, xcoord, ycoord);
      buffer.WriteInteger("%1i",0);
      buffer.NewLine();
    }
    buffer.NewLine();
  }
  // Write buffer
  IO->Write(buffer.Buffer(),buffer.CurrentSize());
  // End and Pause command
  Printf("end\npause -1\n");
  return 0;
} 
