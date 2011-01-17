//DataFile
#include <cstdlib>
#include <cstring>
#include "DataFile.h"
#include "CpptrajStdio.h"

#define DATABUFFERSIZE 128

// CONSTRUCTOR - Arg is datafile name
DataFile::DataFile(char *nameIn) {
  // Default xlabel value
  xlabel=(char*) malloc( 6 * sizeof(char));
  strcpy(xlabel,"Frame");
  noEmptyFrames=false;
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

// Set DataFile debug level
// NOTE: Pass in constructor?
void DataFile::SetDebug(int debugIn) {
  if (debug==debugIn) return;
  debug=debugIn;
  if (debug>0)
    mprintf("DataFile %s DEBUG LEVEL SET TO %i\n",filename,debug);
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
  char buffer[DATABUFFERSIZE];

  // Datafile Header
  outfile->IO->Printf("#%-7s",xlabel);
  for (set=0; set<Nsets; set++) {
    // Skip those empty sets
    if ( SetList[set]->CheckSet() ) continue;
    outfile->IO->Printf(" %12s",SetList[set]->Name());
  }
  outfile->IO->Printf("\n");

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
    outfile->IO->Printf("%8i",frame);
    for (set=0; set<Nsets; set++) {
      // Skip those empty sets
      if ( SetList[set]->CheckSet() ) continue;
      SetList[set]->Write(buffer,frame);
      outfile->IO->Printf("%s",buffer);
    }
    outfile->IO->Printf("\n");
  }
}

/*
 * DataFile::WriteDataInverted()
 * Alternate method of writing out data. Each frame is put into a column, with
 * column 1 containing headers.
 */
void DataFile::WriteDataInverted(PtrajFile *outfile, int maxFrames) {
  int frame,set,empty;
  char buffer[DATABUFFERSIZE];

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
    // Write header as first column
    outfile->IO->Printf("\"%12s\" ",SetList[set]->Name());
    // Write each frame to a column
    for (frame=0; frame<maxFrames; frame++) {
      SetList[set]->Write(buffer,frame);
      outfile->IO->Printf("%s",buffer);
    }
    outfile->IO->Printf("\n");
  }
}
 
/*
 * DataFile::WriteGrace() 
 * Write out sets to file in xmgrace format
 */
void DataFile::WriteGrace(PtrajFile *outfile, int maxFrames) {
  int set,frame;
  char buffer[DATABUFFERSIZE];

  // Grace Header
  outfile->IO->Printf("@with g0\n");
  outfile->IO->Printf("@  xaxis label \"%s\"\n",xlabel);
  outfile->IO->Printf("@  yaxis label \"\"\n");
  outfile->IO->Printf("@  legend 0.2, 0.995\n");
  outfile->IO->Printf("@  legend char size 0.60\n");

  // Loop over sets in data
  for (set=0; set<Nsets; set++) {
    // Skip those empty sets
    if ( SetList[set]->CheckSet() ) continue;
    // Set information
    outfile->IO->Printf("@  s%i legend \"%s\"\n",set,SetList[set]->Name());
    outfile->IO->Printf("@target G0.S%i\n",set);
    outfile->IO->Printf("@type xy\n");
    // Set Data
    for (frame=0; frame<maxFrames; frame++) {
      // If specified, run through every set in the frame and check if empty
      if (noEmptyFrames) {
        if ( SetList[set]->isEmpty(frame) ) continue; 
      }
      outfile->IO->Printf("%8i",frame);
      SetList[set]->Write(buffer,frame);
      outfile->IO->Printf("%s\n",buffer);
    }
  }
}

#undef DATABUFFERSIZE
