//DataFile
#include <cstdlib>
#include <cstring>
#include "DataFile.h"
#include "CpptrajStdio.h"
#include "CharBuffer.h"

/// CONSTRUCTOR
DataFile::DataFile() {
  xlabel=NULL;
  ylabel=NULL;
  noEmptyFrames=false;
  noXcolumn=false;
  filename=NULL;
  SetList=NULL;
  Nsets=0;
  filename=NULL;
  debug=0;
  isInverted=false;
  maxFrames = 0;
  xcol_width = 0;
  x_format = NULL;

  xmin=1;
  xstep=1;
  ymin=1;
  ystep=1;
  useMap=false;
  printLabels=true;
}

/// CONSTRUCTOR - Arg is datafile name
DataFile::DataFile(char *nameIn) {
  // Default xlabel value is Frame
  xlabel=(char*) malloc( 6 * sizeof(char));
  strcpy(xlabel,"Frame");
  xcol_width = 0;
  x_format = NULL;
  // Default ylabel value is blank
  ylabel=(char*) malloc( sizeof(char));
  strcpy(ylabel,"");
  noEmptyFrames=false;
  noXcolumn=false;
  filename=NULL;
  SetList=NULL;
  Nsets=0;
  filename=(char*) malloc( (strlen(nameIn)+1) * sizeof(char));
  strcpy(filename,nameIn);
  debug=0;
  isInverted=false;
  maxFrames = 0;
  xmin=1;
  xstep=1;
  ymin=1;
  ystep=1;
  useMap=false;
  printLabels=true;
}

/// DESTRUCTOR
DataFile::~DataFile() {
  if (filename!=NULL) free(filename);
  // Individual data sets are freed in DataSetList
  if (SetList!=NULL) free(SetList);
  if (xlabel!=NULL) free(xlabel);
  if (ylabel!=NULL) free(ylabel);
  if (x_format!=NULL) delete[] x_format;
}

// DataFile::SetDebug
/// Set DataFile debug level
void DataFile::SetDebug(int debugIn) {
  if (debug==debugIn) return;
  debug=debugIn;
  if (debug>0)
    mprintf("DataFile %s DEBUG LEVEL SET TO %i\n",filename,debug);
}

// DataFile::SetNoXcol
/// Turn printing of frame column off.
void DataFile::SetNoXcol() {
  noXcolumn=true;
}

// DataFile::SetNoEmptyFrames()
/// If called signals that datasets with empty frames should be skipped
void DataFile::SetNoEmptyFrames() {
  noEmptyFrames=true;
}

// DataFile::SetXlabel()
void DataFile::SetXlabel(char *labelIn) {
  if (labelIn==NULL) return;
  xlabel = (char*) realloc ( xlabel, (strlen(labelIn)+1) * sizeof(char) );
  strcpy(xlabel,labelIn);
}

// DataFile::SetYlabel()
void DataFile::SetYlabel(char *labelIn) {
  if (labelIn==NULL) return;
  ylabel = (char*) realloc ( ylabel, (strlen(labelIn)+1) * sizeof(char) );
  strcpy(ylabel,labelIn);
}

// DataFile::SetInverted()
/// If called signals that X and Y values should be flipped (data and grace only)
void DataFile::SetInverted() {
  isInverted=true;
}

// DataFile::SetCoordMinStep()
/// Set the min and step values for the X and Y coords.
void DataFile::SetCoordMinStep(double xminIn, double xstepIn, 
                               double yminIn, double ystepIn) {
  xmin = xminIn;
  xstep = xstepIn;
  ymin = yminIn;
  ystep = ystepIn;
}

// DataFile::SetMap()
/// (Gnuplot only) Turn off the cornerstocolor option
/** gnuplot will be directed to generated a true interpolated contour map.
  */
void DataFile::SetMap() {
  useMap=true;
}

// DataFile::SetNoLabels()
/** Gnuplot only currently. Do not print y axis labels (dataset names). This
  * will automatically be turned off if the # labels > 30.
  */
void DataFile::SetNoLabels() {
  printLabels=false;
}

// DataFile::SetPrecision()
/** Set precision of the specified dataset to width.precision. If '*' specified 
  * set for all datasets in file.
  */
void DataFile::SetPrecision(char *dsetName, int widthIn, int precisionIn) {
  int precision, dset;
  DataSet *Dset = NULL;

  if (dsetName==NULL) {
    mprintf("Error: SetPrecision must be called with dataset name or '*'.\n");
    return;
  }
  if (widthIn<1) {
    mprintf("Error: SetPrecision (%s): Cannot set width < 1.\n",filename);
    return;
  }
  precision=precisionIn;
  if (precisionIn<0) precision=0;
  // If <dsetName>=='*' specified set precision for all data sets
  if (dsetName[0]=='*') {
    mprintf("    Setting width.precision for all sets in %s to %i.%i\n",
            filename,widthIn,precision);
    for (dset=0; dset<Nsets; dset++)
      SetList[dset]->SetPrecision(widthIn,precision);

  // Otherwise find dataset <dsetName> and set precision
  } else {
    mprintf("    Setting width.precision for dataset %s to %i.%i\n",
            dsetName,widthIn,precision);
    for (dset=0; dset<Nsets; dset++) {
      if ( strcmp(SetList[dset]->Name(), dsetName)==0 ) {
        Dset=SetList[dset];
        break;
      }
    }
    if (Dset!=NULL)
      Dset->SetPrecision(widthIn,precision);
    else
      mprintf("Error: Dataset %s not found in datafile %s\n",dsetName,filename);
  }
}

// DataFile::AddSet()
/// Add given set to this datafile
int DataFile::AddSet(DataSet *D) {
  if (D==NULL) return 1;
  SetList = (DataSet**) realloc( SetList, (Nsets+1) * sizeof(DataSet*));
  SetList[Nsets]=D;
  Nsets++;
  return 0;
}

// DataFile::NameIs()
/// Return true if datafile name matches nameIn
bool DataFile::DataFileNameIs(char *nameIn) {
  if (strcmp(nameIn, filename)==0) return true;
  return false;
}

// DataFile::DataSetNames()
/** Print Dataset names to one line. If the number of datasets is greater 
  * than 10 just print the first and last 4 data sets.
  */
void DataFile::DataSetNames() {
  int set;
  if (Nsets>10) {
    for (set=0; set < 4; set++)
      mprintf(" %s",SetList[set]->Name());
    mprintf(" ...");
    for (set=Nsets-4; set < Nsets; set++)
      mprintf(" %s",SetList[set]->Name()); 
    //mprintf("%i datasets.",Nsets);
  } else {
    for (set=0; set<Nsets; set++) {
      if (set>0) mprintf(",");
      mprintf("%s",SetList[set]->Name());
    }
  }
}

// DataFile::Write()
/** Write datasets to file. Check that datasets actually contain data. 
  * Exit if no datasets in this datafile have been used.
  */
void DataFile::Write() {
  CpptrajFile outfile;
  int set,nwrite;
  int maxSetFrames = 0;
  int currentMax = 0;
  int NnewList = 0;
  DataSet **newList;

  // Remove data sets that do not contain data from SetList.
  // NOTE: CheckSet also sets up the format string for the dataset.
  newList = (DataSet**) malloc( Nsets * sizeof(DataSet*) );
  nwrite=0;
  for (set=0; set<Nsets; set++) {
    if ( SetList[set]->CheckSet() ) {
      mprintf("Warning: DataFile %s: Set %s contains no data - skipping.\n",
              filename, SetList[set]->Name());
      //return;
    } else {
      newList[NnewList++] = SetList[set];
      // Determine what the maxium x value for this set is.
      // Should be last value added.
      maxSetFrames = SetList[set]->Xmax();
      if (maxSetFrames > currentMax) currentMax = maxSetFrames;
      nwrite++;
    }
  }
  // Reset SetList
  free(SetList);
  SetList = newList;
  Nsets = NnewList;
  // If all data sets are empty then no need to write
  if (nwrite==0) {
    mprintf("Warning: DataFile %s has no sets containing data - skipping.\n",filename);
    return;
  }
  // Since currentMax is the last frame, increment currentMax by 1 for use in for loops
  maxFrames = currentMax + 1;
  //mprintf("DEBUG: Max frames for %s is %i (maxFrames=%i)\n",filename,currentMax,maxFrames);

  if (outfile.SetupFile(filename,WRITE,UNKNOWN_FORMAT,UNKNOWN_TYPE,debug)) return;
  if (outfile.OpenFile()) return;

  // If format is unknown default to DATAFILE
  if (outfile.fileFormat == UNKNOWN_FORMAT)
    outfile.fileFormat = DATAFILE;

  switch (outfile.fileFormat) {
    case DATAFILE : 
      if (isInverted)
        this->WriteDataInverted(&outfile);
      else 
        this->WriteData(&outfile); 
      break;
    case XMGRACE  : 
      if (isInverted)
        this->WriteGraceInverted(&outfile);
      else 
        this->WriteGrace(&outfile); 
      break;
    case GNUPLOT  :
      if (isInverted)
        mprintf("Warning: Gnuplot format does not support invert; printing standard.\n");
      this->WriteGnuplot(&outfile);
      break;
    default       : mprintf("Error: Datafile %s: Unknown type.\n",filename);
  }

  outfile.CloseFile();
}

// DataFile::WriteData()
/// Write datasets to file. Put each set in its own column. 
void DataFile::WriteData(CpptrajFile *outfile) {
  int set,frame,empty;
  CharBuffer buffer;
  int lineSize = 0;
  bool firstSet = true;
  size_t dataFileSize = 0;
  double xcoord;
  int xcol_precision = 0; // Default x column precision is 0

  // Create format string for X column. If xstep is not 1, set precision to 3
  xcol_width = DigitWidth( maxFrames );
  if (xstep!=1) {
    xcol_precision = 3; 
    // If the width for the x column plus the characters needed for precision
    // (3 + decimal point) would be greater than 8, increment the x column
    // width by 4.
    if ( xcol_width > 4 ) xcol_width+=4; 
  } 
  // Default width for x col is at least 8
  if (xcol_width < 8) xcol_width = 8;
  x_format = SetAlignedDoubleFormatString(xcol_width, xcol_precision);
   
  // Calculate overall size of the datafile. 
  if (!noXcolumn) lineSize = xcol_width;
  for (set=0; set<Nsets; set++)
    lineSize += SetList[set]->Width();
  // Add 1 to lineSize for newline
  lineSize++;
  // Datafile size is line size times maxFrames+1 lines
  dataFileSize = (size_t) maxFrames;
  dataFileSize++;
  dataFileSize *= (size_t) lineSize;
  // NOTE: Since we are manually writing the newlines to this file there is
  //       no need to allocate space for a NULL char. If something like 
  //       sprintf is ever used to write newlines or any data that is at the
  //       end of the file, 1 more byte needs to be allocd for NULL.
  buffer.Allocate( dataFileSize );

  // Write header to buffer
  if (!noXcolumn) {
    // NOTE: Calling WriteStringN with width-1 since it automatically inserts
    //       a space along with the label (before or after according to 
    //       leftAlign). Currently total width written with WriteStringN is
    //       width+1. 
    buffer.WriteStringN(xlabel,xcol_width-1,true); 
  }
  for (set=0; set<Nsets; set++) {
    if (noXcolumn && firstSet) {
      SetList[set]->WriteNameToBuffer(buffer,true);
      firstSet=false;
    } else
      SetList[set]->WriteNameToBuffer(buffer,false);
  }
  buffer.NewLine();

  // Write Data to buffer
  for (frame=0; frame<maxFrames; frame++) {
    // If specified, run through every set in the frame and check if empty
    if (noEmptyFrames) {
      empty=0;
      for (set=0; set<Nsets; set++) {
        if ( SetList[set]->isEmpty(frame) ) {empty++; break;}
      } 
      if (empty!=0) continue;
    }
    // Output Frame
    if (!noXcolumn) {
      xcoord = (xstep * frame) + xmin;
      buffer.WriteDouble(x_format,xcoord);
    }
    for (set=0; set<Nsets; set++) {
      SetList[set]->WriteBuffer(buffer,frame);
    }
    buffer.NewLine();
  }

  // Write buffer to file
  outfile->IO->Write(buffer.Buffer(),sizeof(char),buffer.CurrentSize());
}

// DataFile::WriteDataInverted()
/** Alternate method of writing out data where X and Y values are switched. 
  * Each frame is put into a column, with column 1 containing headers.
  */
void DataFile::WriteDataInverted(CpptrajFile *outfile) {
  int frame,set,empty;
  CharBuffer buffer;
  size_t dataFileSize = 0;

  // Determine max output file size
  // Each line is (set width * (number of frames+1) + newline
  for (set = 0; set < Nsets; set++)
    dataFileSize += ((SetList[set]->Width() * (maxFrames+1)) + 1);
  buffer.Allocate( dataFileSize );

  for (set=0; set < Nsets; set++) {
    // if specified check for empty frames in the set
    if (noEmptyFrames) {
      empty=0;
      for (frame=0; frame<maxFrames; frame++) {
        if (SetList[set]->isEmpty(frame) ) {empty++; break;}
      }
      if (empty!=0) continue;
    }
    // Write dataset name as first column
    SetList[set]->WriteNameToBuffer(buffer,false);
    // Write each frame to subsequent columns
    for (frame=0; frame<maxFrames; frame++) {
      SetList[set]->WriteBuffer(buffer,frame);
    }
    buffer.NewLine();
  }
  outfile->IO->Write(buffer.Buffer(),sizeof(char),buffer.CurrentSize());
}
 
// DataFile::WriteGrace() 
/// Write out sets to file in xmgrace format.
void DataFile::WriteGrace(CpptrajFile *outfile) {
  int set,frame;
  CharBuffer buffer;
  double xcoord;
  size_t dataFileSize;

  // Create format string for X column. Default precision is 3
  xcol_width = DigitWidth( maxFrames );
  // If the width for the x column plus the characters needed for precision
  // (3 + decimal point) would be greater than 8, increment the x column
  // width by 4.
  if ( xcol_width > 4 ) xcol_width+=4;
  // Default width for x col is at least 8
  if (xcol_width < 8) xcol_width = 8;
  x_format = SetAlignedDoubleFormatString(xcol_width, 3);

  // Calculate file size:
  // 1) Initial Header: 91 + xlabel + ylabel
  dataFileSize = (91 + strlen(xlabel) + strlen(ylabel));
  // 2) Set Headers: 37 + 8 + SetName + 8 (fixing integers at size 8)
  //    Check that number of sets will not require > 8 chars
  if (DigitWidth(Nsets) > 8) {
    mprinterr("Internal Error: WriteGrace: # of sets (%i) is too large!\n",Nsets);
    return;
  }
  for (set = 0; set < Nsets; set++) {
    dataFileSize += (53 + strlen( SetList[set]->Name() ));
  // 3) Set Data: (xcol_width + SetWidth + 1) * maxFrames
    dataFileSize += ((xcol_width + SetList[set]->Width() + 1 ) * maxFrames);
  }
  // Allocate buffer
  buffer.Allocate( dataFileSize );

  // Grace header
  buffer.Sprintf("@with g0\n@  xaxis label \"%s\"\n@  yaxis label \"%s\"\n@  legend 0.2, 0.995\n@  legend char size 0.60\n",xlabel,ylabel);

  // Loop over sets in data
  for (set=0; set<Nsets; set++) {
    // Set information
    buffer.Sprintf("@  s%-8i legend \"%s\"\n@target G0.S%-8i\n@type xy\n",
                   set,SetList[set]->Name(),set);
    // Set Data
    for (frame=0; frame<maxFrames; frame++) {
      // If specified, run through every set in the frame and check if empty
      if (noEmptyFrames) {
        if ( SetList[set]->isEmpty(frame) ) continue; 
      }
      xcoord = (xstep * frame) + xmin;
      buffer.WriteDouble(x_format,xcoord);
      SetList[set]->WriteBuffer(buffer,frame);
      buffer.NewLine();
    }
  }
  outfile->IO->Write(buffer.Buffer(),sizeof(char),buffer.CurrentSize());
}

// DataFile::WriteGraceInverted() 
/** Write out sets to file in xmgrace format. Write out data from each
  * frame as 1 set.
  */
void DataFile::WriteGraceInverted(CpptrajFile *outfile) {
  int set,frame,empty;
  size_t dataFileSize;
  CharBuffer buffer;
  double xcoord;

  // Create format string for X column. Default precision is 3
  xcol_width = DigitWidth( maxFrames );
  // If the width for the x column plus the characters needed for precision
  // (3 + decimal point) would be greater than 8, increment the x column
  // width by 4.
  if ( xcol_width > 4 ) xcol_width+=4;
  // Default width for x col is at least 8
  if (xcol_width < 8) xcol_width = 8;
  x_format = SetAlignedDoubleFormatString(xcol_width, 3);

  // Calculate file size
  // 1) Initial header: 91 + xlabel + ylabel
  dataFileSize = (91 + strlen(xlabel) + strlen(ylabel));
  // 2) Set Headers: (22 + 8) * maxFrames  (fixing integers at size 8)
  //    Check that number of frames will not require > 8 chars
  if (DigitWidth(maxFrames) > 8) {
    mprinterr("Internal Error: WriteGraceInverted: # of frames (%i) is too large!\n",maxFrames);
    return;
  }
  dataFileSize += (30 * maxFrames);
  // 3) Set Data: (xcol_width + SetWidth + 4 + SetName) * maxFrames
  for (set = 0; set < Nsets; set++)
    dataFileSize += ((xcol_width + SetList[set]->Width() + 4 + strlen( SetList[set]->Name() ))
                     * maxFrames);
  // Allocate buffer
  buffer.Allocate(dataFileSize);

  // Grace Header
  buffer.Sprintf("@with g0\n@  xaxis label \"%s\"\n@  yaxis label \"%s\"\n@  legend 0.2, 0.995\n@  legend char size 0.60\n",ylabel,xlabel);

  // Loop over frames
  for (frame=0; frame<maxFrames; frame++) {
    // If specified, run through every set in the frame and check if empty
    if (noEmptyFrames) {
      empty=0;
      for (set=0; set<Nsets; set++) {
        if ( SetList[set]->isEmpty(frame) ) {empty++; break;}
      }
      if (empty!=0) continue;
    }

    // Set information
    buffer.Sprintf("@target G0.S%-8i\n@type xy\n",frame);
    // Loop over all Set Data for this frame
    for (set=0; set<Nsets; set++) {
      xcoord = (ystep * set) + ymin;
      buffer.WriteDouble(x_format,xcoord);
      SetList[set]->WriteBuffer(buffer,frame);
      buffer.Sprintf(" \"%s\"",SetList[set]->Name());
      buffer.NewLine();
    }
  }
  outfile->IO->Write(buffer.Buffer(),sizeof(char),buffer.CurrentSize());
}

// DataFile::WriteGnuplot()
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
void DataFile::WriteGnuplot(CpptrajFile *outfile) {
  CharBuffer buffer;
  int set, frame;
  double xcoord, ycoord;
  size_t dataFileSize, xycol, Erow, Eset, dataheader, Set_total;

  // Create format string for X and Y columns. Default precision is 3
  xcol_width = DigitWidth( maxFrames );
  // If the width for the x/y columns plus the characters needed for precision
  // (3 + decimal point) would be greater than 8, increment the column
  // width by 4.
  if ( xcol_width > 4 ) xcol_width+=4;
  // Default width for cols is at least 8
  if (xcol_width < 8) xcol_width = 8;
  x_format = SetAlignedDoubleFormatString(xcol_width, 3);

  // Turn off labels if number of sets is too large since they 
  // become unreadable. Should eventually have some sort of 
  // autotick option.
  if (printLabels && Nsets > 30 ) {
    mprintf("Warning: %s: gnuplot - number of sets %i > 30, turning off Y labels.\n",
            filename,Nsets);
    printLabels=false;
  }

  // Calculate data file size
  // 1) Label, range, plot, and end commands:
  //    87 + xlabel + ylabel + [4*8] + filename + 13
  dataFileSize = 132 + strlen(xlabel) + strlen(ylabel) + strlen(filename);
  // 2) Data labels:
  //    22 + 8 + 8 + 2 + SUM[Nsets]:(4+8+SetName)
  if (printLabels) {
    dataFileSize += 40;
    for (set = 0; set < Nsets; set++)
      dataFileSize += ( 12 + strlen( SetList[set]->Name() ) );
  }
  // 3) Header and data
  //    If not using map gnuplot requires an extra row and column in the X
  //    and Y directions to prevent data from being truncated, so print
  //    an empty row for each set and an empty set after all frames.
  // xycol: Combined size of x and y coords
  // dataheader: Size of the data header (Map command)
  // Erow: size of empty row
  // Eset size of empty set
  xycol = xcol_width + xcol_width;
  if (!useMap) {
    dataheader = 30;
    Erow = xycol + 4;                // 4 == data (%1i) + space + space + newline
    Eset = ((Nsets + 1) * Erow) + 1; // last 1 == newline
  } else {
    dataheader = 13;
    Erow = 0; // No empty row
    Eset = 0; // No empty sets
  }
  //    Compute size of each dataset plus the x and y coords
  Set_total = 0;
  for (set = 0; set < Nsets; set++) 
    Set_total += (xycol + SetList[set]->Width() + 3);
  //    Total size
  dataFileSize += ((maxFrames * ( Set_total + Erow + 1 )) + Eset);
  // Add a little extra in case there is some overflow with setting
  // the range, labels, etc. This is a bit lazy but probably harmless.
  // The real bulk of the data is calculated accurately.
  dataFileSize += 128;

  // Allocate buffer
  buffer.Allocate( dataFileSize );

  // Map command
  if (!useMap)
    buffer.WriteString("%s","set pm3d map corners2color c1\n");
  else
    buffer.WriteString("%s","set pm3d map\n");

  // Y axis Data Labels
  if (printLabels) {
    // NOTE: Add option to turn on grid later?
    //outfile->IO->Printf("set pm3d map hidden3d 100 corners2color c1\n");
    //outfile->IO->Printf("set style line 100 lt 2 lw 0.5\n");
    // Set up Y labels
    buffer.Sprintf("set ytics %8.3f,%8.3f\nset ytics(",ymin,ystep);
    for (set=0; set<Nsets; set++) {
      if (set>0) buffer.WriteString("%s",",");
      ycoord = (ystep * set) + ymin;
      buffer.Sprintf("\"%s\" %8.3f",SetList[set]->Name(),ycoord);
    }
    buffer.Sprintf(")\n");
  }

  // Set axis label and range, write plot command
  // Make Yrange +1 and -1 so entire grid can be seen
  ycoord = (ystep * Nsets) + ymin;
  // Make Xrange +1 and -1 as well
  xcoord = (xstep * maxFrames) + xmin;
  buffer.Sprintf("set xlabel \"%s\"\nset ylabel \"%s\"\nset yrange [%8.3f:%8.3f]\nset xrange [%8.3f:%8.3f]\nsplot \"-\" with pm3d title \"%s\"\n",
                 xlabel,ylabel,ymin-ystep,ycoord+ystep,xmin-xstep,xcoord+xstep,filename);

  // Data
  for (frame=0; frame < maxFrames; frame++) {
    xcoord = (xstep * frame) + xmin;
    for (set=0; set < Nsets; set++) {
      ycoord = (ystep * set) + ymin;
      buffer.WriteDouble(x_format,xcoord);
      buffer.Space();
      buffer.WriteDouble(x_format,ycoord);
      buffer.Space();
      SetList[set]->WriteBuffer(buffer,frame);
      buffer.NewLine();
    }
    if (!useMap) {
      // Print one empty row for gnuplot pm3d without map
      ycoord = (ystep * set) + ymin;
      buffer.WriteDouble(x_format,xcoord);
      buffer.Space();
      buffer.WriteDouble(x_format,ycoord);
      buffer.Space();
      buffer.WriteInteger("%1i",0);
      buffer.NewLine();
    }
    buffer.NewLine();
  }
  if (!useMap) {
    // Print one empty set for gnuplot pm3d without map
    xcoord = (xstep * frame) + xmin;
    for (set=0; set<=Nsets; set++) {
      ycoord = (ystep * set) + ymin;
      buffer.WriteDouble(x_format,xcoord);
      buffer.Space();
      buffer.WriteDouble(x_format,ycoord);
      buffer.Space();
      buffer.WriteInteger("%1i",0);
      buffer.NewLine();
    }
    buffer.NewLine();
  }

  // End and Pause command
  buffer.WriteString("%s","end\npause -1");
  buffer.NewLine();
  // Write buffer
  outfile->IO->Write(buffer.Buffer(),sizeof(char),buffer.CurrentSize());
}

