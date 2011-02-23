//DataFile
#include <cstdlib>
#include <cstring>
#include "DataFile.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR - Arg is datafile name
DataFile::DataFile(char *nameIn) {
  // Default xlabel value is Frame
  xlabel=(char*) malloc( 6 * sizeof(char));
  strcpy(xlabel,"Frame");
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
}

// DESTRUCTOR
DataFile::~DataFile() {
  if (filename!=NULL) free(filename);
  // Individual data sets are freed in DataSetList
  if (SetList!=NULL) free(SetList);
  if (xlabel!=NULL) free(xlabel);
}

/*
 * DataFile::SetDebug
 * Set DataFile debug level
 * NOTE: Pass in constructor?
 */
void DataFile::SetDebug(int debugIn) {
  if (debug==debugIn) return;
  debug=debugIn;
  if (debug>0)
    mprintf("DataFile %s DEBUG LEVEL SET TO %i\n",filename,debug);
}

/*
 * DataFile::SetNoXcol
 * Turn printing of frame column off.
 */
void DataFile::SetNoXcol() {
  noXcolumn=true;
}

/*
 * DataFile::SetXlabel()
 */
void DataFile::SetXlabel(char *labelIn) {
  if (labelIn==NULL) return;
  xlabel = (char*) realloc ( xlabel, (strlen(labelIn)+1) * sizeof(char) );
  strcpy(xlabel,labelIn);
}

/*
 * DataFile::SetInverted()
 */
void DataFile::SetInverted() {
  isInverted=true;
}

/*
 * DataFile::AddSet()
 * Add given set to this datafile
 */
int DataFile::AddSet(DataSet *D) {

  SetList = (DataSet**) realloc( SetList, (Nsets+1) * sizeof(DataSet*));
  SetList[Nsets]=D;
  Nsets++;
  return 0;
}

/* 
 * DataFile::NameIs()
 * Return 1 if datafile name matches nameIn
 */
int DataFile::NameIs(char *nameIn) {
  if (strcmp(nameIn, filename)==0) return 1;
  return 0;
}

/*
 * DataFile::DataSetNames()
 * Print Dataset names to one line
 */
void DataFile::DataSetNames() {
  int set;
  for (set=0; set<Nsets; set++) {
    if (set>0) mprintf(",");
    mprintf("%s",SetList[set]->Name());
  }
}

/*
 * DataFile::Write()
 * Write datasets to file. Check that datasets actually contain data. 
 * Exit if no datasets in this datafile have been used.
 */
void DataFile::Write(bool noEmptyFramesIn) {
  PtrajFile outfile;
  int set,nwrite;
  int maxSetFrames = 0;
  int currentMax = 0;

  noEmptyFrames=noEmptyFramesIn;
  // Check that at least some data sets contain data
  nwrite=0;
  for (set=0; set<Nsets; set++) {
    if ( SetList[set]->CheckSet() ) {
      mprintf("Warning: DataFile %s: Set %s contains no data - skipping.\n",
              filename, SetList[set]->Name());
      //return;
    } else {
      // Determine what the maxium x value for this set is.
      // Should be last value added.
      maxSetFrames = SetList[set]->Xmax();
      if (maxSetFrames > currentMax) currentMax = maxSetFrames;
      nwrite++;
    }
  }
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

/*
 * DataFile::WriteData()
 * Write datasets to file. Put each set in its own column. 
 */
void DataFile::WriteData(PtrajFile *outfile) {
  int set,frame,empty;
  char *ptr,*buffer;
  int lineSize = 0;
  bool firstSet = true;

  buffer=NULL;
  // Datafile Header - Also calculate line size.
  // Add 8 to lineSize for frame if !noXcolumn
  if (!noXcolumn) {
    outfile->IO->Printf("#%-7s",xlabel);
    lineSize += 8;
  }
  for (set=0; set<Nsets; set++) {
    // Skip those empty sets
    if ( SetList[set]->CheckSet() ) continue;
    // Increment lineSize by Width of this set
    lineSize += SetList[set]->Width();
    // Print name of this set
    if (noXcolumn && firstSet) {
      outfile->IO->Printf("#%-12s",SetList[set]->Name());
      firstSet=false;
    } else
      outfile->IO->Printf(" %12s",SetList[set]->Name());
  }
  outfile->IO->Printf("\n");
  // Add 2 to lineSize, 1 for newline, 1 for NULL
  lineSize += 2;
  //mprintf("DEBUG: Calculated lineSize is %i\n",lineSize);
  // Allocate buffer
  buffer = (char*) malloc( lineSize * sizeof(char) );

  // Data
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
    // NOTE: For consistency with Ptraj start at frame 1
    if (!noXcolumn)
      outfile->IO->Printf("%8i",frame + OUTPUTFRAMESHIFT);
    ptr = buffer;
    for (set=0; set<Nsets; set++) {
      // Skip those empty sets
      if ( SetList[set]->CheckSet() ) continue;
      ptr = SetList[set]->Write(ptr,frame);
    }
    // IO->Printf is limited to 1024 chars, use Write instead
    outfile->IO->Write(buffer,sizeof(char),strlen(buffer));
    outfile->IO->Printf("\n");
  }
  // Free buffer
  if (buffer!=NULL) free(buffer);
}

/*
 * DataFile::WriteDataInverted()
 * Alternate method of writing out data where X and Y values are switched. 
 * Each frame is put into a column, with column 1 containing headers.
 */
void DataFile::WriteDataInverted(PtrajFile *outfile) {
  int frame,set,empty;
  int currentLineSize=0;
  int lineSize=0;
  char *buffer,*ptr;

  buffer=NULL;
  for (set=0; set < Nsets; set++) {
    // Skip empty sets
    if (SetList[set]->CheckSet()) continue;
    // if specified check for empty frames in the set
    if (noEmptyFrames) {
      empty=0;
      for (frame=0; frame<maxFrames; frame++) {
        if (SetList[set]->isEmpty(frame) ) {empty++; break;}
      }
      if (empty!=0) continue;
    }
    // Necessary buffer size is (set width * number of frames) + newline + NULL
    lineSize = (SetList[set]->Width() * maxFrames) + 2;
    // Check if lineSize > currentLineSize; if so, reallocate buffer
    if (lineSize > currentLineSize) {
      buffer = (char*) realloc(buffer, lineSize * sizeof(char));
      currentLineSize=lineSize;
    }
    // Write header as first column
    outfile->IO->Printf("\"%12s\" ",SetList[set]->Name());
    // Write each frame to a column
    ptr = buffer;
    for (frame=0; frame<maxFrames; frame++) {
      ptr = SetList[set]->Write(ptr,frame);
    }
    // IO->Printf is limited to 1024 chars, use Write instead
    outfile->IO->Write(buffer,sizeof(char),strlen(buffer));
    outfile->IO->Printf("\n");
  }
  // free buffer
  if (buffer!=NULL) free(buffer);
}
 
/*
 * DataFile::WriteGrace() 
 * Write out sets to file in xmgrace format.
 */
void DataFile::WriteGrace(PtrajFile *outfile) {
  int set,frame;
  int lineSize=0;;
  int currentLineSize=0;
  char *buffer;

  // Grace Header
  outfile->IO->Printf("@with g0\n");
  outfile->IO->Printf("@  xaxis label \"%s\"\n",xlabel);
  outfile->IO->Printf("@  yaxis label \"\"\n");
  outfile->IO->Printf("@  legend 0.2, 0.995\n");
  outfile->IO->Printf("@  legend char size 0.60\n");

  // Loop over sets in data
  buffer=NULL;
  for (set=0; set<Nsets; set++) {
    // Skip those empty sets
    if ( SetList[set]->CheckSet() ) continue;
    // Set information
    outfile->IO->Printf("@  s%i legend \"%s\"\n",set,SetList[set]->Name());
    outfile->IO->Printf("@target G0.S%i\n",set);
    outfile->IO->Printf("@type xy\n");
    // Calc buffer size; reallocate if bigger than current size
    lineSize = SetList[set]->Width() + 2;
    if (lineSize > currentLineSize) {
      buffer = (char*) realloc(buffer, lineSize * sizeof(char));
      currentLineSize = lineSize;
      // Since using IO->Printf, buffer cannot be larger than 1024
      if (lineSize>=1024) {
        mprintf("Error: DataFile::WriteGrace: %s: Data width >= 1024\n",filename);
        if (buffer!=NULL) free(buffer);
        return;
      }
    }
    // Set Data
    for (frame=0; frame<maxFrames; frame++) {
      // If specified, run through every set in the frame and check if empty
      if (noEmptyFrames) {
        if ( SetList[set]->isEmpty(frame) ) continue; 
      }
      SetList[set]->Write(buffer,frame);
      outfile->IO->Printf("%8i%s\n",frame+OUTPUTFRAMESHIFT,buffer);
    }
  }
  if (buffer!=NULL) free(buffer);
}

/*
 * DataFile::WriteGraceInverted() 
 * Write out sets to file in xmgrace format. Write out data from each
 * frame as 1 set.
 */
void DataFile::WriteGraceInverted(PtrajFile *outfile) {
  int set,frame,empty;
  int lineSize=0;;
  int currentLineSize=0;
  char *buffer;

  // Grace Header
  outfile->IO->Printf("@with g0\n");
  outfile->IO->Printf("@  xaxis label \"\"\n");
  outfile->IO->Printf("@  yaxis label \"\"\n");
  outfile->IO->Printf("@  legend 0.2, 0.995\n");
  outfile->IO->Printf("@  legend char size 0.60\n");

  // Calculate maximum expected size for output line.
  for (set=0; set < Nsets; set++) {
    lineSize = SetList[set]->Width() + 2;
    if (lineSize > currentLineSize)
      currentLineSize = lineSize;
  }
  // Since using IO->Printf, buffer cannot be larger than 1024
  if (lineSize>=1024) {
    mprintf("Error: DataFile::WriteGraceInverted: %s: Data width >= 1024\n",filename);
    return;
  }
  buffer = (char*) malloc(lineSize * sizeof(char));

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
    //outfile->IO->Printf("@  s%i legend \"%s\"\n",set,SetList[set]->Name());
    outfile->IO->Printf("@target G0.S%i\n",frame);
    outfile->IO->Printf("@type xy\n");
    // Loop over all Set Data for this frame
    for (set=0; set<Nsets; set++) {
      // Skip those empty sets
      if ( SetList[set]->CheckSet() ) continue;
      SetList[set]->Write(buffer,frame);
      outfile->IO->Printf("%8i%s \"%s\"\n",set+OUTPUTFRAMESHIFT,buffer,SetList[set]->Name());
    }
  }
  if (buffer!=NULL) free(buffer);
}

/*
 * DataFile::WriteGnuplot()
 * Write each frame from all sets in blocks in the following format:
 *   Frame Set   Value
 *   1     0.5   X
 *   1     1.5   X
 *   1     2.5   X
 *
 *   2     0.5   X
 *   2     1.5   X
 *   ...
 * Subtract 0.5 from set number so that grid lines will be centered on Y values
 */
void DataFile::WriteGnuplot(PtrajFile *outfile) {
  char *buffer;
  int set, frame;
  int lineSize=0;
  int currentLineSize=0;
  float Fset;

  // Calculate maximum expected size for output line.
  for (set=0; set < Nsets; set++) {
    lineSize = SetList[set]->Width() + 2;
    if (lineSize > currentLineSize)
      currentLineSize = lineSize;
  }
  // Since using IO->Printf, buffer cannot be larger than 1024
  if (lineSize>=1024) {
    mprintf("Error: DataFile::WriteGnuplot: %s: Data width >= 1024\n",filename);
    return;
  }
  buffer = (char*) malloc(lineSize * sizeof(char));

  // Header
  // Turn off labels if number of sets is too large since they 
  // become unreadable. Should eventually have some sort of 
  // autotick option.
  if (Nsets > 30 ) {
    mprintf("Warning: %s: gnuplot - number of sets %i > 30, turning off Y labels.\n",
            filename,Nsets);
    outfile->IO->Printf("set pm3d map corners2color c1\n");
  } else {
    outfile->IO->Printf("set pm3d map corners2color c1\n");
    // NOTE: Add option to turn on grid later?
    //outfile->IO->Printf("set pm3d map hidden3d 100 corners2color c1\n");
    //outfile->IO->Printf("set style line 100 lt 2 lw 0.5\n");
    // Set up Y labels
    outfile->IO->Printf("set ytics 1,1,%i\nset ytics(",Nsets);
    for (set=0; set<Nsets; set++) {
      if (set>0) outfile->IO->Printf(",");
      outfile->IO->Printf("\"%s\" %i",SetList[set]->Name(),set+1);
    }
    outfile->IO->Printf(")\n");
  }
  outfile->IO->Printf("set xlabel \"%s\"\n",xlabel);
  // Make Yrange +1 and -1 so entire grid can be seen
  outfile->IO->Printf("set yrange [0.0:%i.0]\n",Nsets+1);
  // Make Xrange +1 and -1 as well
  outfile->IO->Printf("set xrange [0.0:%i.0]\n",maxFrames+1);
  // Plot command
  outfile->IO->Printf("splot \"-\" with pm3d title \"%s\"\n",filename);

  // Data
  for (frame=0; frame < maxFrames; frame++) {
    Fset = OUTPUTFRAMESHIFT;
    Fset -= 0.5;
    for (set=0; set < Nsets; set++) {
      SetList[set]->Write(buffer,frame);
      outfile->IO->Printf("%8i %8.1f %s\n",frame+OUTPUTFRAMESHIFT,Fset,buffer);
      Fset++;
    }
    // Print one empty row for gnuplot pm3d
    outfile->IO->Printf("%8i %8.1f %8i\n\n",frame+OUTPUTFRAMESHIFT,Fset,0);
  }

  // End and Pause command
  outfile->IO->Printf("end\npause -1\n");
  // Free buffer
  if (buffer!=NULL) free(buffer);
}

