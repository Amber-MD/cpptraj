#include <cstdio> // sscanf
#include "DataIO_Evecs.h"
#include "CpptrajStdio.h"
#include "BufferedFrame.h"
#include "DataSet_Modes.h"

// CONSTRUCTOR
DataIO_Evecs::DataIO_Evecs() : ibeg_(1), iend_(-1) {
  SetValid( DataSet::MODES );
}

// DataIO_Evecs::ID_DataFormat()
bool DataIO_Evecs::ID_DataFormat(CpptrajFile& infile) {
  if (infile.OpenFile()) return false;
  std::string firstLine = infile.GetLine();
  infile.CloseFile();
  return ( firstLine.compare(0, 18," Eigenvector file:") == 0 );
}

// DataIO_Evecs::ReadHelp()
void DataIO_Evecs::ReadHelp() {
  mprintf("\tibeg <firstmode>: First mode to read in (default 1).\n"
          "\tiend <lastmode>:  Last mode to read in (default 50).\n");
}

// DataIO_Evecs::processReadArgs()
int DataIO_Evecs::processReadArgs(ArgList& argIn) {
  ibeg_ = argIn.getKeyInt("ibeg",1);
  iend_ = argIn.getKeyInt("iend",-1);
  if ((iend_ != -1 && iend_ < 1) || ibeg_ < 1) {
    mprinterr("Error: iend and ibeg must be > 0\n");
    return 1;
  }
  if (iend_ != -1 && iend_ < ibeg_) {
    mprinterr("Error: iend cannot be less than ibeg\n");
    return 1;
  }
  return 0;
}

/** This routine should be updated when new matrix types are added. */
const char* DataIO_Evecs::MatrixOutputString(MetaData::scalarType typeIn) {
  const char* ptr = 0;
  switch (typeIn) {
    case MetaData::DIST :      ptr = "DIST"; break;
    case MetaData::COVAR :     ptr = "COVAR"; break;
    case MetaData::MWCOVAR :   ptr = "MWCOVAR"; break;
    case MetaData::CORREL :    ptr = "CORREL"; break;
    case MetaData::DISTCOVAR : ptr = "DISTCOVAR"; break;
    case MetaData::IDEA :      ptr = "IDEA"; break;
    case MetaData::IREDMAT  :  ptr = "IRED"; break;
    case MetaData::DIHCOVAR :  ptr = "DIHCOVAR"; break;
    default:                  ptr = "UNKNOWN";
  }
  return ptr;
}

// DataIO_Evecs::ReadData()
int DataIO_Evecs::ReadData(FileName const& modesfile,
                           DataSetList& datasetlist, std::string const& dsname)
{
  mprintf("\tReading modes from %s\n", modesfile.full());
  BufferedFrame infile;
  if (infile.OpenRead( modesfile)) return 1;
  // Read title line, convert to arg list
  const char* buffer = 0;
  if ( (buffer = infile.NextLine())==0 ) {
    mprinterr("Error: ReadEvecFile(): error while reading title (%s)\n",infile.Filename().full());
    return 1;
  }
  ArgList title(buffer);
  // Check if reduced
  bool reduced = title.hasKey("Reduced");
  // Determine modes file type
  MetaData::scalarType modesType = MetaData::UNDEFINED;
  for (int matidx = (int)MetaData::DIST;
           matidx != (int)MetaData::UNDEFINED; ++matidx)
  {
    if (title.hasKey( MatrixOutputString((MetaData::scalarType)matidx) ))
    {
      modesType = (MetaData::scalarType)matidx;
      break;
    }
  }
  // For compatibility with quasih and nmode output
  if (modesType == MetaData::UNDEFINED) {
    mprintf("Warning: ReadEvecFile(): Unrecognized type [%s]. Assuming MWCOVAR.\n",
            title.ArgLine());
    modesType = MetaData::MWCOVAR;
  }
  // Allocate MODES dataset. No appending allowed.
  MetaData md(dsname, MetaData::M_MATRIX, modesType);
  DataSet* mds = datasetlist.AddSet( DataSet::MODES, md, "Evecs" );
  if (mds == 0) return 1;
  DataSet_Modes& modesData = static_cast<DataSet_Modes&>( *mds );
  // For newer modesfiles, get # of modes in file.
  int modesInFile = title.getKeyInt("nmodes",-1);
  // Determine number of modes to read
  int modesToRead = 0;
  if (modesInFile == -1) {
    // Older file. If iend not yet set, try default of 50 modes.
    if (iend_ == -1) iend_ = 50;
    modesToRead = iend_ - ibeg_ + 1;
    mprintf("Warning: Older modes file, # of modes not known.\n"
            "Warning: Will try to read at least %i modes.\n", modesToRead);
  } else {
    // Number of modes is known.
    mprintf("\tFile contains %i modes.\n", modesInFile);
    if (iend_ == -1) iend_ = modesInFile;
    if (iend_ > modesInFile) {
      mprintf("Warning: Last mode to read (iend=%i) > modes in file. Setting to %i\n",
              iend_, modesInFile);
      iend_ = modesInFile;
    }
    modesToRead = iend_ - ibeg_ + 1;
  }
  if (modesToRead < 1) {
    mprinterr("Error: No modes to read.\n");
    return 1;
  } else if (modesToRead == 1)
    mprintf("\tAttempting to read mode %i from %s\n", ibeg_, modesfile.full());
  else
    mprintf("\tAttempting to read %i modes (%i to %i) from %s\n", modesToRead,
            ibeg_, iend_, modesfile.full());
  // For newer modesfiles, get width of data elts
  int colwidth = title.getKeyInt("width", -1);
  if (colwidth == -1)
    colwidth = 11; // Default, 10 + 1 space
  modesData.SetupFormat().SetFormatWidthPrecision(colwidth - 1, 5);
  // Read number of elements in avg coords and eigenvectors
  if ( (buffer = infile.NextLine())==0 ) {
    mprinterr("Error: ReadEvecFile(): error while reading number of atoms (%s)\n",
              infile.Filename().full());
    return 1;
  }
  int navgcrd = 0;
  int vecsize = 0;
  int nvals = sscanf(buffer, "%i %i", &navgcrd, &vecsize);
  if (nvals == 0) {
    mprinterr("Error: ReadEvecFile(): sscanf on coords failed (%s)\n",infile.Filename().full());
    return 1;
  } else if (nvals == 1) {
    mprintf("Warning: ReadEvecFile(): No value for eigenvector size found in %s,\n",
            infile.Filename().full());
    mprintf("         assuming it is equal to #average elements (%i)\n",navgcrd);
    vecsize = navgcrd;
  }
  // Allocate FrameBuffer
  int bufsize;
  if (navgcrd > vecsize)
    bufsize = navgcrd;
  else
    bufsize = vecsize;
  infile.SetupFrameBuffer( bufsize, colwidth, 7 );
  // Allocate memory for avg coords and read in
  modesData.AllocateAvgCoords( navgcrd );
  if (navgcrd > 0) {
    infile.ReadFrame( );
    infile.BufferToDouble( modesData.AvgFramePtr(), modesData.NavgCrd() );
    infile.BufferBegin(); // Reset buffer position
  }
  // Allocate buffer memory for eigenvalues and eigenvectors
  double* evalues = new double[ modesToRead ];
  double* evectors = 0;
  if (vecsize > 0)
    evectors = new double[ modesToRead * vecsize ];
  int nmodes = 0;      // Mode currently in memory
  int currentMode = 0; // Modes currently reading in from file
  bool firstRead = true;
  int error_status = 0;
  while ( (buffer = infile.NextLine())!=0 ) { // This should read in ' ****'
    if (buffer[0] != ' ' || buffer[1] != '*' || buffer[2] != '*') {
      mprinterr("Error: ReadEvecFile(): When reading eigenvector %i, expected ' ****',\n",
                currentMode+1);
      mprinterr("Error: got %s [%s]\n", buffer, infile.Filename().full());
      error_status = 1;
      break;
    }
    // Read number and eigenvalue 
    if ( (buffer = infile.NextLine())==0 ) {
      mprinterr("Error: ReadEvecFile(): error while reading number and eigenvalue (%s)\n",
                infile.Filename().full());
      error_status = 2;
      break;
    }
    if (sscanf(buffer, "%*i%lf", evalues + nmodes) != 1) {
      mprinterr("Error: ReadEvecFile(): error while scanning number and eigenvalue (%s)\n",
                infile.Filename().full());
      error_status = 3;
      break;
    }
    if (vecsize > 0) {
      // Read eigenvector
      // Older modesfiles could have vecsize > 0 but no eigenvectors, only 
      // blanks. If this is the case set vecsize to -1 to indicate a blank
      // read is needed after reading eigenvalues.
      double* Vec = evectors + (nmodes * vecsize);
      int vi = 0;
      while (vi < vecsize) {
        buffer = infile.NextLine();
        if (firstRead && (buffer[0] == '\n' || buffer[0] == '\r')) {
          mprintf("Warning: Old modes file with vecsize > 0 but no eigenvectors.\n");
          vecsize = -1;
          delete[] evectors;
          evectors = 0;
          break;
        }
        double tmpval[7];
        int nvals = sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf", tmpval,
                           tmpval+1, tmpval+2, tmpval+3, tmpval+4, tmpval+5, tmpval+6);
        for (int ti = 0; ti < nvals; ++ti)
          Vec[vi++] = tmpval[ti];
      }
      // Check if mode read was between ibeg and iend (which start from 1).
      // If so, increment number of modes.
      if (currentMode+1 >= ibeg_ && currentMode < iend_) ++nmodes;
      if (nmodes == modesToRead) break;
      ++currentMode;
    } else if (vecsize == -1) {
      // Blank read past empty eigenvector
      buffer = infile.NextLine();
    }
    firstRead = false;
  }
  infile.CloseFile();
  if (error_status == 0) {
    if (nmodes != modesToRead)
      mprintf("Warning: Number of read modes is %i, requested %i\n", nmodes, modesToRead);
    error_status = modesData.SetModes( reduced, nmodes, vecsize, evalues, evectors );
  }
  delete[] evalues;
  if (evectors != 0) delete[] evectors;
  return error_status;
}

// DataIO_Evecs::WriteData()
int DataIO_Evecs::WriteData(FileName const& fname, DataSetList const& SetList) {
  if (SetList.empty()) return 1;
  if (SetList.size() > 1)
    mprintf("Warning: Multiple sets not yet supported for Evecs write.\n");
  DataSet_Modes const& modesData = static_cast<DataSet_Modes const&>( *(*(SetList.begin())) );
  BufferedFrame outfile;
  if (outfile.OpenWrite( fname )) {
    mprinterr("Error: Could not open %s for writing.\n", fname.full());
    return 1;
  }
  if (modesData.IsReduced())
    outfile.Printf(" Reduced Eigenvector file: ");
  else
    outfile.Printf(" Eigenvector file: ");
  outfile.Printf("%s", MatrixOutputString(modesData.Meta().ScalarType()));
  // Write out # of modes on title line to not break compat. with older modes files
  outfile.Printf(" nmodes %i", modesData.Nmodes());
  // Write out col width on title line to not break compat. with older modes files
  int colwidth = modesData.Format().ColumnWidth();
  outfile.Printf(" width %i\n", colwidth);
  // First number is # avg coords, second is size of each vector
  outfile.Printf(" %4i %4i\n", modesData.NavgCrd(), modesData.VectorSize());
  // Set up framebuffer, default 7 columns
  int bufsize;
  if (modesData.NavgCrd() > modesData.VectorSize())
    bufsize = modesData.NavgCrd();
  else
    bufsize = modesData.VectorSize();
  outfile.SetupFrameBuffer( bufsize, colwidth, 7 );
  // Print average coords
  outfile.DoubleToBuffer( modesData.AvgFramePtr(), modesData.NavgCrd(),
                          modesData.Format().fmt() );
  outfile.WriteFrame();
  // Eigenvectors and eigenvalues
  for (int mode = 0; mode < modesData.Nmodes(); ++mode) {
    outfile.Printf(" ****\n %4i ", mode+1);
    outfile.Printf(modesData.Format().fmt(), modesData.Eigenvalue(mode));
    outfile.Printf("\n");
    if (modesData.Eigenvectors() != 0) {
      const double* Vec = modesData.Eigenvector(mode);
      outfile.BufferBegin();
      outfile.DoubleToBuffer( Vec, modesData.VectorSize(), modesData.Format().fmt() );
      outfile.WriteFrame();
    }
  }
  outfile.CloseFile();
  return 0;
}
