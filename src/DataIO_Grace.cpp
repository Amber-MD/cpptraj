#include <cstdio> // sscanf
#include <algorithm>
#include "DataIO_Grace.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "DataSet_double.h"
#include "DataSet_string.h"

// TODO: Set dimension labels
// Dont assume anything about set ordering
int DataIO_Grace::ReadData(FileName const& fname, 
                           DataSetList& datasetlist, std::string const& dsname)
{
  ArgList dataline;
  int setnum = 0;
  std::vector<std::string> labels;
  double XY[2];
  const char* linebuffer;
  DataSetList::Darray Xvals;
  DataSetList::DataListType inputSets(1);
  
  // Allocate and set up read buffer
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;

  // Read chunks from file
  linebuffer = buffer.Line(); // First line
  while (linebuffer != 0) {
    if (linebuffer[0] == '@') {
      // Command: create command line without the @
      dataline.SetList(linebuffer+1, " \t");
      if (dataline.empty()) {
        linebuffer = buffer.Line();
        continue;
      }
      if (dataline.CommandIs("target")) {
        // Data to follow. Next line is 'type'
        linebuffer = buffer.Line(); // type line
        if (linebuffer == 0) return 1; // TODO: Better error
        if (linebuffer[0] != '@') return 1; // TODO: Check type
        linebuffer = buffer.Line(); // Should be first line of data.
        DataSet_double* Yvals = new DataSet_double();
        MetaData md(dsname, setnum);
        if ((int)labels.size() > setnum)
          md.SetLegend( labels[setnum] );
        Yvals->SetMeta( md );
        Xvals.clear();
        while (linebuffer != 0 && linebuffer[0] != '@' && 
               linebuffer[0] != '&' && linebuffer[0] != '#')
        {
          if (sscanf(linebuffer, "%lf %lf", XY, XY+1) != 2) break;
          Xvals.push_back(   XY[0] );
          Yvals->AddElement( XY[1] );
          linebuffer = buffer.Line();
        }
        // Should now be positioned 1 line after last data line.
        inputSets[0] = (DataSet*)Yvals;
        if (datasetlist.AddOrAppendSets("", Xvals, inputSets)) return 1;
        ++setnum;
      } else if (dataline[0][0] == 's' || dataline[0][0] == 'S') {
        // Set command
        if (dataline.Nargs() == 3 && dataline[1] == "legend" && !dataline[2].empty())
          labels.push_back(dataline[2]);
        linebuffer = buffer.Line();
      } else
        linebuffer = buffer.Line();
    } else
      linebuffer = buffer.Line();
  }
  buffer.CloseFile();

  // Set DataSet legends if specified
  if (!labels.empty()) {
    if (setnum == (int)labels.size()) {
      mprintf("\tLabels:\n");
      for (unsigned int i = 0; i != labels.size(); i++)
        mprintf("\t\t%s\n", labels[i].c_str());
    } else {
      mprintf("Warning: Number of labels (%zu) does not match number of sets (%i)\n",
              labels.size(), setnum);
    }
  }
  return 0;
}
// -----------------------------------------------------------------------------
void DataIO_Grace::WriteHelp() {
  mprintf("\tinvert: Flip X/Y axes.\n");
  mprintf("\txydy  : Make consecutive sets into XYDY sets.\n");
}

// DataIO_Grace::processWriteArgs()
int DataIO_Grace::processWriteArgs(ArgList &argIn) {
  if (argIn.hasKey("invert")) isInverted_ = true;
  if (argIn.hasKey("xydy")) isXYDY_ = true;
  if (isInverted_ && isXYDY_) {
    mprinterr("Error: 'invert' not compatible with 'xydy'\n");
    return 1;
  }
  return 0;
}

// DataIO_Grace::WriteData()
int DataIO_Grace::WriteData(FileName const& fname, DataSetList const& SetList)
{
  int err = 0;
  // Open output file.
  CpptrajFile file;
  if (file.OpenWrite( fname )) return 1;
  if (isXYDY_)
    err = WriteDataXYDY(file, SetList);
  else if (isInverted_)
    err = WriteDataInverted(file, SetList);
  else
    err = WriteDataNormal(file, SetList);
  file.CloseFile();
  return err;
}

// DataIO_Grace::WriteDataNormal()
int DataIO_Grace::WriteDataNormal(CpptrajFile& file, DataSetList const& Sets) {
  // Hold all 1D data sets.
  if (Sets.empty()) return 1;
  // Determine size of largest DataSet.
  //size_t maxFrames = Sets.DetermineMax();
  // Grace header. Use first data set for labels
  // TODO: DataFile should pass in axis information 
/*  file.Printf(
    "@with g0\n@  xaxis label \"%s\"\n@  yaxis label \"%s\"\n@  legend 0.2, 0.995\n@  legend char size 0.60\n",
    Dim[0].Label().c_str(), Dim[1].Label().c_str()
  );*/
  file.Printf("@with g0\n@  xaxis label \"%s\"\n@  yaxis label \"%s\"\n"
              "@  legend 0.2, 0.995\n@  legend char size 0.60\n",
              Sets[0]->Dim(0).Label().c_str(), "");
  // Determine if we have a single STRING DataSet. Assume these contain labels.
  DataSet_string* labelSet = 0;
  for (DataSetList::const_iterator set = Sets.begin(); set != Sets.end(); ++set)
  {
    if (labelSet == 0) {
      if ( (*set)->Type() == DataSet::STRING )
        labelSet = (DataSet_string*)(*set);
    } else {
      if ( (*set)->Type() == DataSet::STRING ) {
        labelSet = 0;
        break;
      }
    }
  }
  if (labelSet != 0)
    mprintf("\tUsing string set '%s' for data point labels\n", labelSet->legend());
  // Loop over DataSets
  unsigned int setnum = 0;
  DataSet::SizeArray frame(1);
  for (DataSetList::const_iterator set = Sets.begin(); set != Sets.end(); ++set, ++setnum)
  {
    // Skip the label set if defined.
    if (*set == labelSet) continue;
    size_t maxFrames = (*set)->Size();
    // Set information
    file.Printf("@  s%u legend \"%s\"\n@target G0.S%u\n@type xy\n",
                   setnum, (*set)->legend(), setnum );
    // Setup set X coord format.
    TextFormat xfmt(XcolFmt());
    if (XcolPrecSet())
      xfmt = TextFormat( XcolFmt(), XcolWidth(), XcolPrec() );
    else
      xfmt.SetCoordFormat( maxFrames, (*set)->Dim(0).Min(), (*set)->Dim(0).Step(), 8, 3 );
    // Write Data for set
    for (frame[0] = 0; frame[0] < maxFrames; frame[0]++) {
      file.Printf(xfmt.fmt(), (*set)->Coord(0, frame[0]));
      (*set)->WriteBuffer(file, frame);
      if (labelSet != 0)
        file.Printf(" \"%s\"", (*labelSet)[frame[0]].c_str());
      file.Printf("\n");
    }
  }
  return 0;
}

int DataIO_Grace::WriteDataXYDY(CpptrajFile& file, DataSetList const& Sets) {
  // Hold all 1D data sets.
  if (Sets.empty()) return 1;
  // For this mode need even number of sets
  unsigned int maxSets = Sets.size();
  if ((maxSets % 2) != 0) {
    maxSets = maxSets - 1;
    mprintf("Warning: XYDY output requires even number of sets.\n");
    if (maxSets < 1) return 1;
    mprintf("Warning: Only using the first %zu sets.\n", maxSets);
  }
  // Grace header. Use first data set for labels
  // TODO: DataFile should pass in axis information 
  file.Printf("@with g0\n@  xaxis label \"%s\"\n@  yaxis label \"%s\"\n"
              "@  legend 0.2, 0.995\n@  legend char size 0.60\n",
              Sets[0]->Dim(0).Label().c_str(), "");
  // Loop over DataSets
  unsigned int setnum = 0;
  DataSet::SizeArray frame(1);
  for (unsigned int sidx = 0; sidx < maxSets; sidx += 2, ++setnum)
  {
    DataSet* ds1 = Sets[sidx];
    DataSet* ds2 = Sets[sidx+1];
    if (ds1->Size() != ds2->Size())
      mprintf("Warning: Sets %s and %s have different sizes.\n", ds1->legend(), ds2->legend());
    size_t maxFrames = std::min(ds1->Size(), ds2->Size());
    // Set information
    file.Printf("@  s%u legend \"%s\"\n@target G0.S%u\n@type xydy\n",
                   setnum, ds1->legend(), setnum );
    // Setup set X coord format.
    TextFormat xfmt(XcolFmt());
    if (XcolPrecSet())
      xfmt = TextFormat( XcolFmt(), XcolWidth(), XcolPrec() );
    else
      xfmt.SetCoordFormat( maxFrames, ds1->Dim(0).Min(), ds1->Dim(0).Step(), 8, 3 );
    // Write Data for sets
    for (frame[0] = 0; frame[0] < maxFrames; frame[0]++) {
      file.Printf(xfmt.fmt(), ds1->Coord(0, frame[0]));
      ds1->WriteBuffer(file, frame);
      ds2->WriteBuffer(file, frame);
      file.Printf("\n");
    }
  }
  return 0;
}

// DataIO_Grace::WriteDataInverted()
int DataIO_Grace::WriteDataInverted(CpptrajFile& file, DataSetList const& Sets)
{
  // Hold all 1D data sets.
  if (Sets.empty()) return 1;
  // Determine size of largest DataSet.
  size_t maxFrames = DetermineMax( Sets );
  // Grace header. Use first DataSet for axis labels.
  // TODO: DataFile should pass in axis info.
/*  file.Printf(
    "@with g0\n@  xaxis label \"%s\"\n@  yaxis label \"%s\"\n@  legend 0.2, 0.995\n@  legend char size 0.60\n",
    Dim[1].Label().c_str(), Dim[0].Label().c_str() 
  );*/
  file.Printf("@with g0\n@  xaxis label \"%s\"\n@  yaxis label \"%s\"\n"
              "@  legend 0.2, 0.995\n@  legend char size 0.60\n",
              "", Sets[0]->Dim(0).Label().c_str());
  // Setup set X coord format. 
  Dimension Xdim(0.0, 1.0);
  TextFormat xfmt(XcolFmt());
  if (XcolPrecSet())
    xfmt = TextFormat( XcolFmt(), XcolWidth(), XcolPrec() );
  else
    xfmt.SetCoordFormat( Sets.size(), Xdim.Min(), Xdim.Step(), 8, 3 );
  // Loop over frames
  DataSet::SizeArray frame(1);
  for (frame[0] = 0; frame[0] < maxFrames; frame[0]++) {
    // Set information
    file.Printf("@target G0.S%zu\n@type xy\n",frame[0]);
    // Loop over all Set Data for this frame
    unsigned int setnum = 0;
    for (DataSetList::const_iterator set = Sets.begin(); set != Sets.end(); ++set, ++setnum) {
      file.Printf(xfmt.fmt(), Xdim.Coord( setnum ));
      (*set)->WriteBuffer(file, frame);
      file.Printf(" \"%s\"\n", (*set)->legend());
    }
  }
  return 0;
}
