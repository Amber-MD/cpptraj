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
void DataFile::Write(int maxFrames, bool noEmptyFramesIn) {
  PtrajFile outfile;
  int set,nwrite;

  noEmptyFrames=noEmptyFramesIn;
  // Check that at least some data sets contain data
  nwrite=0;
  for (set=0; set<Nsets; set++) {
    if ( SetList[set]->CheckSet() ) {
      mprintf("Warning: DataFile %s: Set %s contains no data - skipping.\n",
              filename, SetList[set]->Name());
      //return;
    } else {
      nwrite++;
    }
  }
  if (nwrite==0) {
    mprintf("Warning: DataFile %s has no sets containing data - skipping.\n",filename);
    return;
  }

  if (outfile.SetupFile(filename,WRITE,UNKNOWN_FORMAT,UNKNOWN_TYPE,debug)) return;
  if (outfile.OpenFile()) return;

  // If format is unknown default to DATAFILE
  if (outfile.fileFormat == UNKNOWN_FORMAT)
    outfile.fileFormat = DATAFILE;

  switch (outfile.fileFormat) {
    case DATAFILE   : 
      if (isInverted)
        this->WriteDataInverted(&outfile,maxFrames);
      else 
        this->WriteData(&outfile, maxFrames); 
      break;
    case XMGRACE    : this->WriteGrace(&outfile, maxFrames); break;
    default      : mprintf("Error: Datafile %s: Unknown type.\n",filename);
  }

  outfile.CloseFile();
}

/*
 * DataFile::WriteData()
 * Write datasets to file. Put each set in its own column. 
 */
void DataFile::WriteData(PtrajFile *outfile, int maxFrames) {
  int set,frame,empty;
  char *buffer;
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
    for (set=0; set<Nsets; set++) {
      // Skip those empty sets
      if ( SetList[set]->CheckSet() ) continue;
      SetList[set]->Write(buffer,frame);
      outfile->IO->Printf("%s",buffer);
    }
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
void DataFile::WriteDataInverted(PtrajFile *outfile, int maxFrames) {
  int frame,set,empty;
  int currentLineSize=0;
  int lineSize=0;
  char *buffer;

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
    for (frame=0; frame<maxFrames; frame++) {
      SetList[set]->Write(buffer,frame);
      outfile->IO->Printf("%s",buffer);
    }
    outfile->IO->Printf("\n");
  }
  // free buffer
  if (buffer!=NULL) free(buffer);
}
 
/*
 * DataFile::WriteGrace() 
 * Write out sets to file in xmgrace format.
 */
void DataFile::WriteGrace(PtrajFile *outfile, int maxFrames) {
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
    }
    // Set Data
    for (frame=0; frame<maxFrames; frame++) {
      // If specified, run through every set in the frame and check if empty
      if (noEmptyFrames) {
        if ( SetList[set]->isEmpty(frame) ) continue; 
      }
      outfile->IO->Printf("%8i",frame + OUTPUTFRAMESHIFT);
      SetList[set]->Write(buffer,frame);
      outfile->IO->Printf("%s\n",buffer);
    }
  }
  if (buffer!=NULL) free(buffer);
}

