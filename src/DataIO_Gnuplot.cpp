#include <cstdio> // sscanf
#include "DataIO_Gnuplot.h"
#include "CpptrajStdio.h"
#include "Array1D.h"
#include "DataSet_2D.h"

// CONSTRUCTOR
DataIO_Gnuplot::DataIO_Gnuplot() :
  DataIO(true, true, false), // Valid for 1D and 2D 
  pm3d_(C2C),
  printLabels_(true),
  useMap_(false),
  jpegout_(false),
  binary_(false),
  writeHeader_(true)
{}

bool DataIO_Gnuplot::ID_DataFormat(CpptrajFile& infile) {
  bool isGnuplot = false;
/*  if (!infile.OpenFile()) {
    ArgList line( infile.NextLine(), " " );
    std::string setarg = line.GetStringKey("set");
    if (setarg == "size" || setarg == "pm3d") isGnuplot = true;
    infile.CloseFile();
  }*/
  return isGnuplot;
}

int DataIO_Gnuplot::ReadData(FileName const& fname,
                             DataSetList& DSL, std::string const& dsname)
{
  // Simple test to determine if this is binary data or not.
  CpptrajFile infile;
  if (infile.OpenRead(fname)) return 1;
  const char* ptr = infile.NextLine();
  if (ptr == 0) return 1;
  char firstChar = ptr[0];
  infile.CloseFile();
  if (firstChar == 0)
    return ReadBinaryData(fname, DSL, dsname, "X", "Y");
  else
    return ReadAsciiHeader(fname, DSL, dsname);
}

int DataIO_Gnuplot::ReadBinaryData(FileName const& fname,
                                  DataSetList& DSL, std::string const& dsname,
                                  std::string const& xlabel, std::string const& ylabel)
{
  mprintf("\tGnuplot data appears to be in binary format.\n");
  CpptrajFile bIn;
  if (bIn.OpenRead(fname)) return 1;
  typedef std::vector<double> Darray;
  Darray Xvals, Yvals;
  Darray matrix_Rmajor;

  // Read number of columns
  float fval;
  bIn.Read( &fval, sizeof(float) );
  int ncols = (int)fval;
  // Read X values and convert to double.
  std::vector<float> Vals( ncols );
  Xvals.reserve( ncols );
  bIn.Read( &Vals[0], ncols*sizeof(float) );
  for (std::vector<float>::const_iterator it = Vals.begin(); it != Vals.end(); ++it)
    Xvals.push_back( (double)*it );
  // Keep reading rows until no more found.
  int nread = 1;
  while (nread > 0) {
    // Read Y value
    if (bIn.Read( &fval, sizeof(float) ) != sizeof(float)) break;
    Yvals.push_back( (double)fval );
    bIn.Read( &Vals[0], ncols*sizeof(float) );
    for (std::vector<float>::const_iterator it = Vals.begin(); it != Vals.end(); ++it)
      matrix_Rmajor.push_back( (double)*it );
  }
  bIn.CloseFile();
  mprintf("\t%zu rows, %i cols (%zu), %zu vals\n", Yvals.size(), ncols, Xvals.size(), matrix_Rmajor.size());

  DataSet* ds = DetermineMatrixType( matrix_Rmajor, Yvals.size(), ncols, DSL, dsname );
  if (ds == 0) return 1;
  Dimension Xdim, Ydim;
  if (!Xdim.SetDimension(Xvals, xlabel))
    mprintf("Warning: X dimension is NOT monotonic.\n");
  if (!Ydim.SetDimension(Yvals, ylabel))
    mprintf("Warning: Y dimension is NOT monotonic.\n");
  ds->SetDim(Dimension::X, Xdim);
  ds->SetDim(Dimension::Y, Ydim);

  return 0;
}

/** NOTE: This assumes data is in the file and is an splot "heat map" */
int DataIO_Gnuplot::ReadAsciiHeader(FileName const& fname,
                                    DataSetList& DSL, std::string const& dsname)
{
  mprintf("\tReading Gnuplot header.\n");
  BufferedLine infile;
  if (infile.OpenFileRead( fname )) return 1;
  // Get past all of the 'set' lines
  std::string xlabel, ylabel;
  const char* ptr = infile.Line();
  while (ptr != 0 && ptr[0] == 's' && ptr[1] == 'e' && ptr[2] == 't') {
    if ( (ptr[4] == 'x' || ptr[4] == 'y') && ptr[5] == 'l' ) {
      ArgList line( ptr, " " );
      if (ptr[4] == 'x')
        xlabel = line.GetStringKey("xlabel");
      else
        ylabel = line.GetStringKey("ylabel");
    }
    ptr = infile.Line();
  }
  if (ptr == 0) {
    mprinterr("Error: No data detected in Gnuplot file.\n");
    return 1;
  }
  // Search for 'splot' command.
  while (ptr != 0 && ptr[0] != 's' && ptr[1] != 'p' && ptr[2] != 'l' &&
                     ptr[3] != 'o' && ptr[4] != 't')
    ptr = infile.Line();
  if (ptr == 0) {
    mprinterr("Error: 'splot' not found in '%s'. CPPTRAJ currently only reads\n"
              "Error:   CPPTRAJ-style Gnuplot files.\n", fname.full());
    return 1;
  }
  // Determine if data is in this file or another file. If it is another file,
  // assume the data is binary.
  ArgList splot_line(ptr, " ");
  std::string splot_name = splot_line.GetStringKey("splot");
  int err = 0;
  if ( splot_name == "-" )
    err = ReadAsciiData(infile, DSL, dsname, xlabel, ylabel);
  else
    err = ReadBinaryData(splot_name, DSL, dsname, xlabel, ylabel);

  infile.CloseFile();
  return err;
}

int DataIO_Gnuplot::ReadAsciiData(BufferedLine& infile,
                                  DataSetList& DSL, std::string const& dsname,
                                  std::string const& xlabel, std::string const& ylabel)
{
  mprintf("\tGnuplot data appears to be in ASCII format.\n");
  // Allocate full matrix. Assume column-major order.
  typedef std::vector<double> Darray;
  Darray Xvals;
  Darray Yvals;
  Darray matrix_Cmajor;
  int ncols = 0;
  int nrows = -1;
  const char* ptr = infile.Line(); // First line of data.
  while (ptr != 0 && ptr[0] != 'e' && ptr[1] != 'n' && ptr[2] != 'd') {
    int row = 0;
    while (ptr != 0 && ptr[0] != '\0') {
      double xval, yval, val;
      sscanf(ptr, "%lf %lf %lf", &xval, &yval, &val);
      //mprintf("DEBUG: '%i' %f %f %f (%i, %i)\n", (int)ptr[0], xval, yval, val, ncols,row);
      if (row   == 0) Xvals.push_back( xval );
      if (ncols == 0) Yvals.push_back( yval );
      matrix_Cmajor.push_back( val );
      row++;
      ptr = infile.Line();
    }
    if (nrows == -1)
      nrows = row;
    else if (nrows != row) {
      mprinterr("Error: Number of rows has changed from %i to %i (line %i)\n",
                nrows, row, infile.LineNumber());
      return 1;
    }
    ptr = infile.Line();
    ncols++;
  }
  mprintf("\t%i rows (%zu), %i cols (%zu), %zu values\n", nrows, Yvals.size(), ncols, Xvals.size(), matrix_Cmajor.size());
  // DEBUG
  if (debug_ > 0) {
    mprintf("Xvals:");
    for (Darray::const_iterator it = Xvals.begin(); it != Xvals.end(); ++it)
      mprintf(" %f", *it);
    mprintf("\nYvals:");
    for (Darray::const_iterator it = Yvals.begin(); it != Yvals.end(); ++it)
      mprintf(" %f", *it);
    mprintf("\nVals:\n");
  }
  // Need to convert column major to row major for DetermineMatrixType
  Darray matrix_Rmajor;
  matrix_Rmajor.reserve( matrix_Cmajor.size() );
  for (int r = 0; r != nrows; r++) {
    int idx = r;
    for (int c = 0; c != ncols; c++, idx += nrows) {
      if (debug_ > 0) mprintf(" %g", matrix_Cmajor[idx]);
      matrix_Rmajor.push_back( matrix_Cmajor[idx] );
    }
    if (debug_ > 0) mprintf("\n");
  }
  matrix_Cmajor.clear();

  DataSet* ds = DetermineMatrixType( matrix_Rmajor, nrows, ncols, DSL, dsname );
  if (ds == 0) return 1;
  Dimension Xdim, Ydim;
  if (!Xdim.SetDimension(Xvals, xlabel))
    mprintf("Warning: X dimension is NOT monotonic.\n");
  if (!Ydim.SetDimension(Yvals, ylabel))
    mprintf("Warning: Y dimension is NOT monotonic.\n");
  ds->SetDim(Dimension::X, Xdim);
  ds->SetDim(Dimension::Y, Ydim);

  return 0;
}

// -----------------------------------------------------------------------------
DataIO_Gnuplot::LabelArray DataIO_Gnuplot::LabelArg( std::string const& labelarg) 
{
  return ArgList( labelarg, ",").List();
}

void DataIO_Gnuplot::WriteHelp() {
  mprintf("\tnolabels       : Do not print axis labels.\n"
          "\tusemap         : pm3d output with 1 extra empty row/col (may improve look).\n"
          "\tpm3d           : Normal pm3d map output.\n"
          "\tnopm3d         : Turn off pm3d\n"
          "\tjpeg           : Plot will write to a JPEG file when used with gnuplot.\n"
          "\ttitle          : Plot title. Default is file name.\n"
//          "\tbinary:   Use binary output\n"
          "\tnoheader       : Do not format plot; data output only.\n"
          "\ttitle <title>  : Set plot title (default file base name).\n"
          "\tpalette <arg>  : Change gnuplot pm3d palette to <arg>:\n"
          "\t          'rgb'   - Red, yellow, green, cyan, blue, magenta, red.\n"
          "\t          'kbvyw' - Black, blue, violet, yellow, white.\n"
          "\t          'bgyr'  - Blue, green, yellow, red.\n"
          "\t          'gray'  - Grayscale.\n"
          "\txlabels <list> : Set x axis labels with comma-separated list, e.g.\n"
          "\t                     'xlabels X1,X2,X3'\n"
          "\tylabels <list> : Set y axis labels.\n"
          "\tzlabels <list> : Set z axis labels.\n");
}

// DataIO_Gnuplot::processWriteArgs()
int DataIO_Gnuplot::processWriteArgs(ArgList &argIn) {
  if (argIn.hasKey("nolabels")) printLabels_ = false;
  if (argIn.hasKey("usemap")) pm3d_ = MAP;
  if (argIn.hasKey("pm3d")) pm3d_ = ON;
  if (argIn.hasKey("nopm3d")) pm3d_ = OFF;
  if (argIn.hasKey("jpeg")) jpegout_ = true;
  if (argIn.hasKey("binary")) binary_ = true;
  if (argIn.hasKey("noheader")) writeHeader_ = false;
  title_ = argIn.GetStringKey("title");
  if (!writeHeader_ && jpegout_) {
    mprintf("Warning: jpeg output not supported with 'noheader' option.\n");
    jpegout_ = false;
  }
  title_ = argIn.GetStringKey("title");
  palette_ = argIn.GetStringKey("palette");
  if (!palette_.empty()) {
    if (pm3d_ == OFF) {
      mprintf("Warning: 'palette' not used when 'nopm3d' specified.\n");
      palette_.clear();
    } else if (palette_ == "rgb")
      palette_.assign("set palette model HSV\nset palette rgb 3,2,2\n");
    else if (palette_ == "kbvyw")
      palette_.assign("set palette rgb 30,31,32\n");
    else if (palette_ == "bgyr")
      palette_.assign("set palette rgb 33,13,10\n");
    else if (palette_ == "gray")
      palette_.assign("set pal gray\n");
    else {
      mprintf("Warning: Unrecognized palette '%s'; ignoring.\n", palette_.c_str());
      palette_.clear();
    }
  }
  // Label arguments
  if (argIn.Contains("xlabels"))
    Xlabels_ = LabelArg( argIn.GetStringKey( "xlabels" ) );
  if (argIn.Contains("ylabels"))
    Ylabels_ = LabelArg( argIn.GetStringKey( "ylabels" ) );
  if (argIn.Contains("zlabels"))
    Zlabels_ = LabelArg( argIn.GetStringKey( "zlabels" ) );
  if (pm3d_ == MAP) useMap_ = true;
  return 0;
}

// DataIO_Gnuplot::Pm3d()
/** Set up command for gnuplot pm3d */
std::string DataIO_Gnuplot::Pm3d(size_t maxFrames) {
  // PM3D command
  std::string pm3d_cmd = "with pm3d";
  switch (pm3d_) {
    case C2C: 
      if (maxFrames == 1)
        file_.Printf("set pm3d map corners2color c3\n");
      else 
        file_.Printf("set pm3d map corners2color c1\n"); 
      break;
    case MAP: file_.Printf("set pm3d map\n"); break;
    case ON : file_.Printf("set pm3d\n"); break;
    case OFF: pm3d_cmd.clear(); break;
  }
  if (!pm3d_cmd.empty() && !palette_.empty())
    file_.Printf("%s", palette_.c_str());
  return pm3d_cmd;
}

// DataIO_Gnuplot::WriteRangeAndHeader()
/** Write gnuplot range, plot labels, and plot command. */
void DataIO_Gnuplot::WriteRangeAndHeader(Dimension const& Xdim, size_t Xmax,
                                         Dimension const& Ydim, size_t Ymax,
                                         std::string const& pm3dstr)
{
  const char* binaryFlag = " ";
  if (binary_) binaryFlag = " binary ";
  file_.Printf("set xlabel \"%s\"\nset ylabel \"%s\"\n", 
               Xdim.Label().c_str(), Ydim.Label().c_str());
  file_.Printf("set yrange [%8.3f:%8.3f]\nset xrange [%8.3f:%8.3f]\n", 
         Ydim.Coord(0) - Ydim.Step(), Ydim.Coord(Ymax + 1),
         Xdim.Coord(0) - Xdim.Step(), Xdim.Coord(Xmax + 1));
  const char* tout;
  if (!title_.empty())
    tout = title_.c_str();
  else
    tout = file_.Filename().base();
  file_.Printf("splot \"%s\"%s%s title \"%s\"\n", data_fname_.full(), binaryFlag, pm3dstr.c_str(),
               tout);
}

// DataIO_Gnuplot::Finish()
void DataIO_Gnuplot::Finish() {
  if (!jpegout_ && writeHeader_) {
    if (binary_)
      file_.Printf("pause -1\n");
    else
      file_.Printf("end\npause -1\n");
  }
}

// DataIO_Gnuplot::JpegOut()
/** Write commands to direct gnuplot to print directly to JPEG. */
void DataIO_Gnuplot::JpegOut(size_t xsize, size_t ysize) {
  if (jpegout_) {
    std::string sizearg = "1024,768";
    // For now, if xsize == ysize make square, otherwise make rectangle.
    if (xsize == ysize)
      sizearg = "768,768";
    // Create jpg filename
    std::string jpegname = file_.Filename().Full() + ".jpg";
    file_.Printf("set terminal jpeg size %s\nset output \"%s\"\n",
                  sizearg.c_str(), jpegname.c_str());
  } else {
    // If not writing jpeg and xsize == ysize, make output square
    if (xsize == ysize)
      file_.Printf("set size square\n");
  }
}

int DataIO_Gnuplot::WriteData(FileName const& fname, DataSetList const& SetList) 
{
  // Determine 1D or 2D.
  if (SetList.empty()) return 0;
  int err = 0;
  // Open output file
  if (writeHeader_ || !binary_) {
    if (file_.OpenWrite( fname )) return 1;
  }
  if (binary_) {
    // If writing binary, determine file name.
    data_fname_ = fname;
    if (writeHeader_)
      data_fname_.Append( ".data" );
  } else {
    // Not binary. Data will be included in plot.
    data_fname_.SetFileName_NoExpansion("-");
  }
  // One dimension
  if (SetList[0]->Ndim() == 1) {
    err = WriteSets1D( SetList );
  } else if (SetList[0]->Ndim() == 2) {
    // Warn about writing multiple sets
    if (SetList.size() > 1)
      mprintf("Warning: %s: Writing multiple 2D sets in GNUplot format may"
              " result in unexpected behavior\n", fname.full());
    for (DataSetList::const_iterator set = SetList.begin(); set != SetList.end(); ++set)
      err += WriteSet2D( *(*set) );
  } else
    err = 1;
  if (file_.IsOpen()) file_.CloseFile();
  return err;
}

const char* DataIO_Gnuplot::BasicPalette[]= {
  "#000000", // Black, 0
  "#0000FF", // Blue,  1
  "#00FF00", // Green, N/2
  "#FF0000", // Red,   N
};

// DataIO_Gnuplot::WriteDefinedPalette()
/** Write out a defined palette to the gnuplot file. */
void DataIO_Gnuplot::WriteDefinedPalette(int ncolors) {
  float mincolor = -0.5;
  float maxcolor = (float)ncolors - 0.5;
  file_.Printf("set cbrange [%8.3f:%8.3f]\nset cbtics %8.3f %8.3f 1.0\n",
               mincolor, maxcolor, mincolor + 0.5, maxcolor - 0.5);
  file_.Printf("set palette maxcolors %i\n", ncolors);
  // NOTE: Giving gnuplot too many colors can mess up the palette 
  //       interpolation, leading to unwanted colors being inserted.
  //       Instead, just define a few "hint" colors; the zero color,
  //       then low/middle/high.
  const char** CurrentPalette = BasicPalette;
  file_.Printf("set palette defined (");
  file_.Printf("0 \"%s\",", CurrentPalette[0]);
  file_.Printf("1 \"%s\",", CurrentPalette[1]);
  if (ncolors > 3)
    file_.Printf("%i \"%s\",", (ncolors / 2), CurrentPalette[2]);
  file_.Printf("%i \"%s\")\n", (ncolors - 1), CurrentPalette[3]);
}

void DataIO_Gnuplot::WriteLabels(LabelArray const& labels, Dimension const& dim, char xyz)
{
  unsigned int offset = 1;
  // If the number of labels is greater than max, write fewer labels.
  const unsigned int max_labels = 20;
  if ( labels.size() > max_labels ) 
    offset = labels.size() / max_labels;
  file_.Printf("set %ctics %8.3f,%8.3f\nset %ctics(", xyz, dim.Coord(0),
               dim.Step() * (double)offset, xyz);
  for (unsigned int i = 0; i < labels.size(); i += offset) {
    if (i>0) file_.Printf(",");
    file_.Printf("\"%s\" %8.3f", labels[i].c_str(), dim.Coord(i));
  }
  file_.Printf(")\n");
}

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
int DataIO_Gnuplot::WriteSets1D(DataSetList const& Sets)
{
  // FIXME: Check that dimension of each set matches.
  if (Sets.empty()) return 1;
  CheckXDimension( Sets );
  // Determine size of largest DataSet.
  size_t maxFrames = DetermineMax( Sets );
  // Use X dimension of set 0 for all set dimensions.
  DataSet* Xdata = Sets[0];
  Dimension const& Xdim = static_cast<Dimension const&>( Xdata->Dim(0) ); 
  Dimension Ydim( 1.0, 1.0 );
  TextFormat x_format(XcolFmt()), y_format(XcolFmt());
  if (XcolPrecSet()) {
    x_format = TextFormat(XcolFmt(), XcolWidth(), XcolPrec());
    y_format = x_format;
  } else {
    x_format.SetCoordFormat( maxFrames,   Xdim.Min(), Xdim.Step(), 8, 3 );
    y_format.SetCoordFormat( Sets.size(), Ydim.Min(), Ydim.Step(), 8, 3 );
  }
  std::string xyfmt = x_format.Fmt() + " " + y_format.Fmt() + " ";

  // Turn off labels if number of sets is too large since they 
  // become unreadable. Should eventually have some sort of 
  // autotick option.
  if (writeHeader_) {
    // Check for JPEG output
    JpegOut( maxFrames, Sets.size() );

    // PM3D command
    std::string pm3d_cmd = Pm3d(maxFrames);

    // Y axis Data Labels
    if (printLabels_) {
      // NOTE: Add option to turn on grid later?
      //outfile->file_.Printf("set pm3d map hidden3d 100 corners2color c1\n");
      //outfile->file_.Printf("set style line 100 lt 2 lw 0.5\n");
      // Set up Y labels
      file_.Printf("set ytics %8.3f,%8.3f\nset ytics(", Ydim.Min(), Ydim.Step());
      std::string label_fmt = "\"%s\" " + y_format.Fmt();
      for (unsigned int iy = 0; iy != Sets.size(); iy++) {
        if (iy>0) file_.Printf(",");
        file_.Printf(label_fmt.c_str(), Sets[iy]->legend(), Ydim.Coord(iy));
      }
      file_.Printf(")\n");
      // Set up Z labels
      if (!Zlabels_.empty()) {
        WriteDefinedPalette(Zlabels_.size());
        file_.Printf("set cbtics(");
        int iz = 0;
        for (std::vector<std::string>::iterator label = Zlabels_.begin();
                                                label != Zlabels_.end(); ++label)
        {
          if (iz > 0) file_.Printf(",");
          file_.Printf("\"%s\" %8.3f", (*label).c_str(), (float)iz++);
        }
        file_.Printf(")\n");
      }
    }
    // Set axis label and range, write plot command
    // Make Yrange +1 and -1 so entire grid can be seen
    WriteRangeAndHeader(Xdim, maxFrames, Ydim, Sets.size(), pm3d_cmd);
  }

  // Data
  if (binary_) {
    // ----- BINARY FORMAT -----------------------
    // TODO !useMap_
    // Only support writing scalar 1D sets.
    std::vector<DataSet_1D*> Bsets;
    for (DataSetList::const_iterator set = Sets.begin(); set != Sets.end(); ++set) {
      if ( (*set)->Group() == DataSet::SCALAR_1D )
        Bsets.push_back( (DataSet_1D*)*set );
      else
        mprintf("Warning: Set '%s' is not 1D scalar; cannot be written in Gnuplot binary format.\n",
                (*set)->legend());
    }
    CpptrajFile bOut;
    if (bOut.OpenWrite( data_fname_ )) return 1;
    mprintf("\tWriting binary gnuplot data to '%s'\n", data_fname_.full());
    // Write number of frames (columns)
    float fval = (float)maxFrames;
    bOut.Write( &fval, sizeof(float) );
    // Convert X values (i.e. columns) to floats and write.
    std::vector<float> Vals( maxFrames );
    for (unsigned int i = 0; i != maxFrames; i++)
      Vals[i] = (float)Xdata->Coord(0, i);
    bOut.Write( &Vals[0], maxFrames*sizeof(float) );
    // For each set (row), write Y value and all data.
    for (unsigned int setnum = 0; setnum != Bsets.size(); setnum++)
    {
      fval = (float)Ydim.Coord(setnum);
      bOut.Write( &fval, sizeof(float) );
      for (unsigned int i = 0; i != maxFrames; i++)
        Vals[i] = Bsets[setnum]->Dval(i);
      bOut.Write( &Vals[0], maxFrames*sizeof(float) );
    }
    bOut.CloseFile();
  } else {
    // ----- ASCII FORMAT ------------------------
    DataSet::SizeArray frame(1, 0);
    for (frame[0] = 0; frame[0] < maxFrames; frame[0]++) {
      double xcoord = Xdata->Coord(0, frame[0]);
      for (size_t setnum = 0; setnum < Sets.size(); ++setnum) {
        file_.Printf( xyfmt.c_str(), xcoord, Ydim.Coord(setnum) );
        Sets[setnum]->WriteBuffer( file_, frame );
        file_.Printf("\n");
      }
      if (!useMap_) {
        // Print one empty row for gnuplot pm3d without map
        file_.Printf( xyfmt.c_str(), xcoord, Ydim.Coord(Sets.size()) );
        file_.Printf("0\n");
      }
      file_.Printf("\n");
    }
    if (!useMap_) {
      // Print one empty set for gnuplot pm3d without map
      double xcoord = Xdata->Coord(0, maxFrames);
      for (size_t blankset=0; blankset <= Sets.size(); blankset++) {
        file_.Printf( xyfmt.c_str(), xcoord, Ydim.Coord(blankset) );
        file_.Printf("0\n");
      }
      file_.Printf("\n");
    }
  }
  Finish();
  return 0;
}

// DataIO_Gnuplot::WriteSet2D()
int DataIO_Gnuplot::WriteSet2D( DataSet const& setIn ) {
  if (setIn.Ndim() != 2) {
    mprinterr("Internal Error: DataSet %s in DataFile %s has %zu dimensions, expected 2.\n",
              setIn.legend(), file_.Filename().full(), setIn.Ndim());
    return 1;
  }
  DataSet_2D const& set = static_cast<DataSet_2D const&>( setIn );

  Dimension const& Xdim = setIn.Dim(0);
  Dimension const& Ydim = setIn.Dim(1);
  if (writeHeader_) {
    // Check for JPEG output
    JpegOut( set.Ncols(), set.Nrows() );

    // PM3D command
    std::string pm3d_cmd = Pm3d(set.Size());

    // Axes Data Labels
    if (printLabels_) {
      // Set up X and Y labels
      if (!Ylabels_.empty()) {
        if ( Ylabels_.size() != set.Nrows() )
          mprintf("Warning: # of Ylabels (%zu) does not match Y dimension (%u)\n",
                  Ylabels_.size(), set.Nrows());
        WriteLabels(Ylabels_, Ydim, 'y');
      }
      if (!Xlabels_.empty()) {
        if ( Xlabels_.size() != set.Ncols() )
          mprintf("Warning: # of Xlabels (%zu) does not match X dimension (%u)\n",
                  Xlabels_.size(), set.Ncols());
        WriteLabels(Xlabels_, Xdim, 'x'); 
      }
    }
    // Set axis label and range, write plot command
    // Make Yrange +1 and -1 so entire grid can be seen
    WriteRangeAndHeader(Xdim, set.Ncols(), Ydim, set.Nrows(), pm3d_cmd);
  }
  if (binary_ ) {
    // ----- BINARY FORMAT -----------------------
    CpptrajFile bOut;
    if (bOut.OpenWrite( data_fname_ )) return 1;
    mprintf("\tWriting binary gnuplot data to '%s'\n", data_fname_.full());
    // Write number of frames (columns)
    float fval = (float)set.Ncols();
    bOut.Write( &fval, sizeof(float) );
    // Convert X values (i.e. columns) to floats and write.
    std::vector<float> Vals( set.Ncols() );
    for (unsigned int i = 0; i != set.Ncols(); i++)
      Vals[i] = (float)set.Coord(0, i);
    bOut.Write( &Vals[0], set.Ncols()*sizeof(float) );
    // For each set (row), write Y value and all data.
    for (unsigned int j = 0; j != set.Nrows(); j++)
    {
      fval = (float)set.Coord(1, j);
      bOut.Write( &fval, sizeof(float) );
      for (unsigned int i = 0; i != set.Ncols(); i++)
        Vals[i] = (float)set.GetElement(i, j);
      bOut.Write( &Vals[0], set.Ncols()*sizeof(float) );
    }
    bOut.CloseFile();
  } else {
    // ----- ASCII FORMAT ------------------------
    // Setup XY coord format
    TextFormat x_fmt(XcolFmt()), y_fmt(XcolFmt());
    if (XcolPrecSet()) {
      x_fmt = TextFormat(XcolFmt(), XcolWidth(), XcolPrec());
      y_fmt = x_fmt;
    } else {
      x_fmt.SetCoordFormat( set.Ncols(), Xdim.Min(), Xdim.Step(), 8, 3 );
      y_fmt.SetCoordFormat( set.Nrows(), Ydim.Min(), Ydim.Step(), 8, 3 );
    }
    std::string xyfmt = x_fmt.Fmt() + " " + y_fmt.Fmt(); // FIXME No trailing space for bkwds compat

    DataSet::SizeArray positions(2, 0);
    for (positions[0] = 0; positions[0] < set.Ncols(); ++positions[0]) {
      double xcoord = set.Coord(0, positions[0]);
      for (positions[1] = 0; positions[1] < set.Nrows(); ++positions[1]) {
        file_.Printf( xyfmt.c_str(), xcoord, set.Coord(1, positions[1]) );
        set.WriteBuffer( file_, positions );
        file_.Printf("\n");
      }
      if (!useMap_) {
        // Print one empty row for gnuplot pm3d without map
        file_.Printf( xyfmt.c_str(), xcoord, set.Coord(1, set.Nrows()) );
        file_.Printf(" 0\n");
      }
      file_.Printf("\n");
    }
    if (!useMap_) {
      // Print one empty set for gnuplot pm3d without map
      double xcoord = set.Coord(0, set.Ncols());
      for (size_t blankset=0; blankset <= set.Nrows(); ++blankset) {
        file_.Printf( xyfmt.c_str(), xcoord, set.Coord(1, blankset) );
        file_.Printf(" 0\n");
      }
      file_.Printf("\n");
    }
  }
  Finish();

  return 0;
}
