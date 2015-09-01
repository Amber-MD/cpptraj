#include <cstdio> // sscanf
#include "DataIO_Grace.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"

// TODO: Set dimension labels
// Dont assume anything about set ordering
int DataIO_Grace::ReadData(FileName const& fname, 
                           DataSetList& datasetlist, std::string const& dsname)
{
  ArgList dataline;
  int setnum = 0;
  std::vector<std::string> labels;
  double XY[0];
  const char* linebuffer;
  DataSetList::DataListType inputSets;
  
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
        DataSet* ds = 0;
        MetaData md( dsname, setnum );
        ds = datasetlist.CheckForSet( md );
        if (ds == 0) {
          // Create new data set.
          ds = datasetlist.AddSet( DataSet::XYMESH, md );
          if (ds == 0) return 1;
        } else {
          // Appending to existing. Only XYMESH allowed.
          if (ds->Type() != DataSet::XYMESH) {
            mprinterr("Error: Append currently requires existing set '%s' to be XYMESH\n",
                      ds->legend());
            return 1;
          }
        }
        int Ndata = 0;
        while (linebuffer != 0 && linebuffer[0] != '@' && 
               linebuffer[0] != '&' && linebuffer[0] != '#')
        {
          if (sscanf(linebuffer, "%lf %lf", XY, XY+1) != 2) break;
          ds->Add(Ndata++, XY);
          linebuffer = buffer.Line();
        }
        // Should now be positioned 1 line after last data line.
        inputSets.push_back( ds );
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
    if (inputSets.size() == labels.size()) {
      mprintf("\tLabels:\n");
      for (unsigned int i = 0; i != labels.size(); i++) {
        mprintf("\t\t%s\n", labels[i].c_str());
        inputSets[i]->SetLegend( labels[i] );
      }
    } else {
      mprintf("Warning: Number of labels (%zu) does not match number of sets (%zu)\n",
              labels.size(), inputSets.size());
    }
  }
  return 0;
}
// -----------------------------------------------------------------------------
void DataIO_Grace::WriteHelp() {
  mprintf("\tinvert: Flip X/Y axes.\n");
}

// DataIO_Grace::processWriteArgs()
int DataIO_Grace::processWriteArgs(ArgList &argIn) {
  isInverted_ = argIn.hasKey("invert");
  return 0;
}

// DataIO_Grace::WriteData()
int DataIO_Grace::WriteData(FileName const& fname, DataSetList const& SetList)
{
  int err = 0;
  // Open output file.
  CpptrajFile file;
  if (file.OpenWrite( fname )) return 1;
  if (isInverted_)
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
  // Loop over DataSets
  unsigned int setnum = 0;
  DataSet::SizeArray frame(1);
  for (DataSetList::const_iterator set = Sets.begin(); set != Sets.end(); ++set, ++setnum) {
    size_t maxFrames = (*set)->Size();
    // Set information
    file.Printf("@  s%u legend \"%s\"\n@target G0.S%u\n@type xy\n",
                   setnum, (*set)->legend(), setnum );
    // Setup set X coord format.
    TextFormat xfmt;
    xfmt.SetCoordFormat( maxFrames, (*set)->Dim(0).Min(), (*set)->Dim(0).Step(), 8, 3 );
    // For backwards compat. introduce a leading space
    (*set)->SetupFormat().SetFormatAlign( TextFormat::LEADING_SPACE );
    // Write Data for set
    for (frame[0] = 0; frame[0] < maxFrames; frame[0]++) {
      file.Printf(xfmt.fmt(), (*set)->Coord(0, frame[0]));
      (*set)->WriteBuffer(file, frame);
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
  Dimension Xdim;
  Xdim.SetStep( 1 );
  TextFormat xfmt;
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
