#include "DataIO_Gnuplot.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataIO_Gnuplot::DataIO_Gnuplot() {
  y_label_ = "";
  ymin_ = 1;
  ystep_ = 1;
  pm3d_ = C2C;
  printLabels_ = true;
  useMap_ = false;
  jpegout_ = false;
  binary_ = false;
}

void DataIO_Gnuplot::LabelArg( std::vector<std::string>& labels, std::string const& labelarg) 
{
  if (!labelarg.empty()) {
    ArgList commasep( labelarg, "," );
    for (int i = 0; i < commasep.Nargs(); ++i)
      labels.push_back( commasep[i] );
  }
}

// DataIO_Gnuplot::processWriteArgs()
int DataIO_Gnuplot::processWriteArgs(ArgList &argIn) {
  ArgList::ConstArg ylabel = argIn.getKeyString("ylabel");
  if (ylabel!=NULL) y_label_.assign(ylabel);
  ymin_ = argIn.getKeyDouble("ymin",ymin_);
  ystep_ = argIn.getKeyDouble("ystep",ystep_);

  if (argIn.hasKey("nolabels")) printLabels_ = false;
  if (argIn.hasKey("usemap")) pm3d_ = MAP;
  if (argIn.hasKey("pm3d")) pm3d_ = ON;
  if (argIn.hasKey("nopm3d")) pm3d_ = OFF;
  if (argIn.hasKey("jpeg")) jpegout_ = true;
  if (argIn.hasKey("binary")) binary_ = true;

  // Label arguments
  LabelArg( Xlabels_, argIn.GetStringKey( "xlabels" ) );
  LabelArg( Ylabels_, argIn.GetStringKey( "ylabels" ) );

  if (pm3d_ == MAP) useMap_ = true;
  return 0;
}

// DataIO_Gnuplot::Pm3d()
/** Set up command for gnuplot pm3d */
std::string DataIO_Gnuplot::Pm3d() {
  // PM3D command
  std::string pm3d_cmd = "with pm3d";
  switch (pm3d_) {
    case C2C: 
      if (maxFrames_ == 1)
        Printf("set pm3d map corners2color c3\n");
      else 
        Printf("set pm3d map corners2color c1\n"); 
      break;
    case MAP: Printf("set pm3d map\n"); break;
    case ON : Printf("set pm3d\n"); break;
    case OFF: pm3d_cmd.clear(); break;
  }
  return pm3d_cmd;
}

// DataIO_Gnuplot::WriteRangeAndHeader()
/** Write gnuplot range, plot labels, and plot command. */
void DataIO_Gnuplot::WriteRangeAndHeader(double xcoord, double ycoord, 
                                         std::string const& pm3dstr)
{
  Printf("set xlabel \"%s\"\nset ylabel \"%s\"\n", x_label_.c_str(), y_label_.c_str());
  Printf("set yrange [%8.3f:%8.3f]\nset xrange [%8.3f:%8.3f]\n", 
         ymin_ - ystep_, ycoord + ystep_,
         xmin_ - xstep_, xcoord + xstep_);
  Printf("splot \"-\" %s title \"%s\"\n", pm3dstr.c_str(), BaseFileStr());
}

// DataIO_Gnuplot::Finish()
void DataIO_Gnuplot::Finish() {
  if (!jpegout_)
    Printf("end\npause -1\n");
}

// DataIO_Gnuplot::JpegOut()
/** Write commands to direct gnuplot to print directly to JPEG. */
void DataIO_Gnuplot::JpegOut(int xsize, int ysize) {
  if (jpegout_) {
    std::string sizearg = "1024,768";
    // For now, if xsize == ysize make square, otherwise make rectangle.
    if (xsize == ysize)
      sizearg = "768,768";
    // Create jpg filename
    std::string jpegname = FullFileName() + ".jpg";
    Printf("set terminal jpeg size %s\nset output \"%s\"\n",
                  sizearg.c_str(), jpegname.c_str());
  }
}

int DataIO_Gnuplot::WriteData(DataSetList &SetList) {
  //mprintf("BINARY IS %i\n", (int)binary_);
  if (binary_)
    return WriteDataBinary( SetList );
  else
    return WriteDataAscii( SetList );
}

/** Format:
  *   <N+1>  <y0>   <y1>   <y2>  ...  <yN>
  *    <x0> <z0,0> <z0,1> <z0,2> ... <z0,N>
  *    <x1> <z1,0> <z1,1> <z1,2> ... <z1,N>
  *     :      :      :      :   ...    :
  */
int DataIO_Gnuplot::WriteDataBinary(DataSetList &SetList) {
  DataSetList::const_iterator set;

  int Ymax = SetList.size();
  if (!useMap_)
    ++Ymax;
  float fvar = (float)Ymax;
  mprintf("Ymax = %f\n",fvar);
  IO->Write( &fvar, sizeof(float) );
  for (int setnum = 0; setnum < Ymax; ++setnum) {
    double ycoord = (ystep_ * (double)setnum) + ymin_;
    fvar = (float)ycoord;
    IO->Write( &fvar, sizeof(float) );
  }
  // Data
  for (int frame = 0; frame < maxFrames_; frame++) {
    double xcoord = (xstep_ * (double)frame) + xmin_;
    fvar = (float)xcoord;
    IO->Write( &fvar, sizeof(float) );
    for (set=SetList.begin(); set !=SetList.end(); set++) {
      fvar = (float)(*set)->Dval( frame );
      IO->Write( &fvar, sizeof(float) );
    }
    if (!useMap_) {
      // Print one empty row for gnuplot pm3d without map
      fvar = 0;
      IO->Write( &fvar, sizeof(float) );
    }
  }
  if (!useMap_) {
    // Print one empty set for gnuplot pm3d without map
    double xcoord = (xstep_ * (double)maxFrames_) + xmin_;
    fvar = (float)xcoord;
    IO->Write( &fvar, sizeof(float) );
    fvar = 0;
    for (int blankset=0; blankset < Ymax; blankset++)
      IO->Write( &fvar, sizeof(float) ); 
  }

  return 0;
}

// DataIO_Gnuplot::WriteData()
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
int DataIO_Gnuplot::WriteDataAscii(DataSetList &SetList) {
  DataSetList::const_iterator set;
  double xcoord, ycoord;

  // Create format string for X and Y columns. Default precision is 3
  SetupXcolumn();
  std::string xy_format_string = x_format_ + " " + x_format_ + " ";
  const char *xy_format = xy_format_string.c_str();

  // Turn off labels if number of sets is too large since they 
  // become unreadable. Should eventually have some sort of 
  // autotick option.
/*  if (printLabels_ && SetList.size() > 30 ) {
    mprintf("Warning: %s: gnuplot: number of sets (%i) > 30, turning off Y labels.\n",
            BaseName(), SetList.size());
    printLabels_ = false;
  }*/

  // Check for JPEG output
  JpegOut( maxFrames_, (int)SetList.size() );

  // PM3D command
  std::string pm3d_cmd = Pm3d();

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
      Printf("\"%s\" %8.3f",(*set)->Legend().c_str(),ycoord);
      ++setnum; 
    }
    Printf(")\n");
  }

  // Set axis label and range, write plot command
  // Make Yrange +1 and -1 so entire grid can be seen
  ycoord = (ystep_ * SetList.size()) + ymin_;
  // Make Xrange +1 and -1 as well
  xcoord = (xstep_ * maxFrames_) + xmin_;

  WriteRangeAndHeader(xcoord, ycoord, pm3d_cmd);

  // Data
  int frame = 0;
  for (; frame < maxFrames_; frame++) {
    // If not printing empty frames, make sure that every set has data
    // at this frame.
    if (!printEmptyFrames_) {
      bool emptyFrames = false;
      for (set = SetList.begin(); set != SetList.end(); set++) {
        if ( (*set)->FrameIsEmpty(frame) ) {
          emptyFrames = true;
          break;
        }
      }
      if (emptyFrames) continue;
    }
    xcoord = (xstep_ * frame) + xmin_;
    int setnum = 0;
    for (set=SetList.begin(); set !=SetList.end(); set++) {
      ycoord = (ystep_ * setnum) + ymin_;
      Printf( xy_format, xcoord, ycoord );
      (*set)->WriteBuffer( *this, frame );
      Printf("\n");
      ++setnum;
    }
    if (!useMap_) {
      // Print one empty row for gnuplot pm3d without map
      ycoord = (ystep_ * setnum) + ymin_;
      Printf(xy_format,xcoord,ycoord);
      Printf("0\n");
    }
    Printf("\n");
  }
  if (!useMap_) {
    // Print one empty set for gnuplot pm3d without map
    xcoord = (xstep_ * frame) + xmin_;
    for (int blankset=0; blankset <= (int)SetList.size(); blankset++) {
      ycoord = (ystep_ * blankset) + ymin_;
      Printf(xy_format,xcoord,ycoord);
      Printf("0\n");
    }
    Printf("\n");
  }
  // End and Pause command
  Finish();
  return 0;
}

// DataIO_Gnuplot::WriteData2D()
int DataIO_Gnuplot::WriteData2D( DataSet& set ) {
  std::vector<int> dimensions;

  // Get dimensions
  set.GetDimensions(dimensions);
  if (dimensions.size() != 2) {
    mprinterr("Internal Error: DataSet %s in DataFile %s has %zu dimensions, expected 2.\n",
              set.c_str(), FullFileStr(), dimensions.size());
    return 1;
  } 

  // Check for JPEG output
  JpegOut( dimensions[0], dimensions[1] );

  // PM3D command
  std::string pm3d_cmd = Pm3d();

  // Axes Data Labels
  if (printLabels_) {
    // Set up X and Y labels
    if (!Ylabels_.empty()) {
      if ( (int)Ylabels_.size() != dimensions[1])
        mprintf("Warning: # of Ylabels (%zu) does not match Y dimension (%i)\n",
                Ylabels_.size(), dimensions[1]);
      Printf("set ytics %8.3f,%8.3f\nset ytics(",ymin_,ystep_);
      for (int iy = 0; iy < (int)Ylabels_.size(); ++iy) {
        if (iy>0) Printf(",");
        double ycoord = (ystep_ * (double)iy) + ymin_;
        Printf("\"%s\" %8.3f", Ylabels_[iy].c_str(), ycoord);
      }
      Printf(")\n");
    }
    if (!Xlabels_.empty()) {
      if ( (int)Xlabels_.size() != dimensions[0])
        mprintf("Warning: # of Xlabels (%zu) does not match X dimension (%i)\n",
                Xlabels_.size(), dimensions[0]); 
      Printf("set xtics %8.3f,%8.3f\nset xtics(",xmin_,xstep_);
      for (int ix = 0; ix < (int)Xlabels_.size(); ++ix) {
        if (ix>0) Printf(",");
        double xcoord = (xstep_ * (double)ix) + xmin_;
        Printf("\"%s\" %8.3f", Xlabels_[ix].c_str(), xcoord);
      }
      Printf(")\n");
    }
  }

  // Set axis label and range, write plot command
  // Make Yrange +1 and -1 so entire grid can be seen
  double ycoord = (ystep_ * (double)dimensions[1]) + ymin_;
  // Make Xrange +1 and -1 as well
  double xcoord = (xstep_ * (double)dimensions[0]) + xmin_;
  WriteRangeAndHeader(xcoord, ycoord, pm3d_cmd);

  for (int ix = 0; ix < dimensions[0]; ++ix) {
    double xcoord = (xstep_ * (double)ix) + xmin_;
    for (int iy = 0; iy < dimensions[1]; ++iy) {
      double ycoord = (ystep_ * (double)iy) + ymin_;
      Printf("%8.3f %8.3f", xcoord, ycoord);
      set.Write2D( *this, ix, iy );
      Printf("\n");
    }
    if (!useMap_) {
      // Print one empty row for gnuplot pm3d without map
      ycoord = (ystep_ * dimensions[1]) + ymin_;
      Printf("%8.3f %8.3f 0\n", xcoord, ycoord);
    }
    Printf("\n");
  }
  if (!useMap_) {
    // Print one empty set for gnuplot pm3d without map
    xcoord = (xstep_ * dimensions[0]) + xmin_;
    for (int blankset=0; blankset <= dimensions[1]; blankset++) {
      ycoord = (ystep_ * blankset) + ymin_;
      Printf("%8.3f %8.3f 0\n", xcoord, ycoord);
    }
    Printf("\n");
  }
  // End and Pause command
  Finish();

  return 0;
}
