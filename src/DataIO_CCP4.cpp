#include <cmath> //sqrt
#include "DataIO_CCP4.h"
#include "CpptrajStdio.h"
#include "ByteRoutines.h"
#include "DataSet_GridFlt.h"

const size_t DataIO_CCP4::wSize = 4;

bool DataIO_CCP4::MapCharsValid(const unsigned char* MAP) {
  return (MAP[0] == 'M' && MAP[1] == 'A' && MAP[2] == 'P' && MAP[3] == ' ');
}

// DataIO_CCP4::ID_DataFormat()
bool DataIO_CCP4::ID_DataFormat( CpptrajFile& infile ) {
  bool isCCP4 = false;
  if (!infile.OpenFile()) {
    unsigned char MAP[4];
    MAP[0] = ' '; MAP[1] = ' '; MAP[2] = ' '; MAP[3] = ' ';
    if (infile.Seek(52 * wSize) == 0) {
      infile.Read( MAP, wSize );
      isCCP4 = MapCharsValid( MAP );
    }
    infile.CloseFile();
  }
  return isCCP4;
}

// DataIO_CCP4::ReadData()
/** Header is 256 4-byte words. Integer unless otherwise noted. First 56 words are:
  *    0-2: columns, rows, sections (fastest changing to slowest)
  *      3: mode: 0 = envelope stored as signed bytes (from -128 lowest to 127 highest)
  *               1 = Image     stored as Integer*2
  *               2 = Image     stored as Reals
  *               3 = Transform stored as Complex Integer*2
  *               4 = Transform stored as Complex Reals
  *               5 == 0
  *    4-6: Column, row, and section offsets
  *    7-9: Intervals along X, Y, Z
  *  10-15: float; 3x cell lengths (Ang) and 3x cell angles (deg)
  *  16-18: Map of which axes correspond to cols, rows, sections (1,2,3 = x,y,z)
  *  19-21: float; Min, max, and mean density
  *  22-24: Space group, bytes used for storing symm ops, flag for skew transform
  *         If skew flag != 0, skew transformation is from standard orthogonal
  *         coordinate frame (as used for atoms) to orthogonal map frame, as:
  *             Xo(map) = S * (Xo(atoms) - t)
  *  25-33: Skew matrix 'S' (in order S11, S12, S13, S21 etc)
  *  34-36: Skew translation 't'
  *  37-51: For future use and can be skipped.
  *     52: char; 'MAP '
  *     53: char; machine stamp for determining endianness
  *     54: float; RMS deviation of map from mean
  *     55: Number of labels
  */
int DataIO_CCP4::ReadData(FileName const& fname,
                            DataSetList& datasetlist, std::string const& dsname)
{
  CpptrajFile infile;
  if (infile.OpenRead( fname )) return 1;
  // Read first 56 words of the header into a buffer.
  headerbyte buffer;
  if (infile.Read(buffer.i, 224*sizeof(unsigned char)) < 1) {
    mprinterr("Error: Could not buffer CCP4 header.\n");
    return 1;
  }
  if (debug_ > 0)
    mprintf("DEBUG: MAP= '%c %c %c %c'  MACHST= '%x %x %x %x'\n",
            buffer.c[208], buffer.c[209], buffer.c[210], buffer.c[211],
            buffer.c[212], buffer.c[213], buffer.c[214], buffer.c[215]);
  // SANITY CHECK
  if (!MapCharsValid(buffer.c + 208)) {
    mprinterr("Error: CCP4 file missing 'MAP ' string at word 53\n");
    return 1;
  }
  // Check endianess
  bool isBigEndian = (buffer.c[212] == 0x11 && buffer.c[213] == 0x11 &&
                      buffer.c[214] == 0x00 && buffer.c[215] == 0x00);
  if (!isBigEndian) {
    if (debug_ > 0) mprintf("DEBUG: Little endian.\n");
    // SANITY CHECK
    if ( !(buffer.c[212] == 0x44 && buffer.c[213] == 0x41 &&
           buffer.c[214] == 0x00 && buffer.c[215] == 0x00) )
      mprintf("Warning: Invalid machine stamp: %x %x %x %x : assuming little endian.\n",
              buffer.c[212], buffer.c[213], buffer.c[214], buffer.c[215]);
  } else {
    if (debug_ > 0) mprintf("DEBUG: Big endian.\n");
    // Perform endian swapping on header if necessary
    endian_swap(buffer.i, 56);
  }

  // Print DEBUG info
  if (debug_ > 0) {
    mprintf("DEBUG: Columns=%i  Rows=%i  Sections=%i\n", buffer.i[0], buffer.i[1], buffer.i[2]);
    mprintf("DEBUG: Mode=%i\n", buffer.i[3]);
    mprintf("DEBUG: Offsets: C=%i  R=%i  S=%i\n", buffer.i[4], buffer.i[5], buffer.i[6]);
    mprintf("DEBUG: NXYZ={ %i %i %i }\n", buffer.i[7], buffer.i[8], buffer.i[9]);
    mprintf("DEBUG: Box XYZ={ %f %f %f }  ABG={ %f %f %f }\n",
            buffer.f[10], buffer.f[11], buffer.f[12],
            buffer.f[13], buffer.f[14], buffer.f[15]);
    mprintf("DEBUG: Map: ColAxis=%i  RowAxis=%i  SecAxis=%i\n",
            buffer.i[16], buffer.i[17], buffer.i[18]);
    mprintf("DEBUG: SpaceGroup#=%i  SymmOpBytes=%i  SkewFlag=%i\n",
            buffer.i[22], buffer.i[23], buffer.i[24]);
    const int* MSKEW = buffer.i + 25;
    mprintf("DEBUG: Skew matrix: %i %i %i\n"
            "                    %i %i %i\n"
            "                    %i %i %i\n", MSKEW[0], MSKEW[1], MSKEW[2], MSKEW[3],
            MSKEW[4], MSKEW[5], MSKEW[6], MSKEW[7], MSKEW[8]);
    const int* TSKEW = buffer.i + 34;
    mprintf("DEBUG: Skew translation: %i %i %i\n", TSKEW[0], TSKEW[1], TSKEW[2]);
    mprintf("DEBUG: Nlabels=%i\n", buffer.i[55]);
  }
  // Check input data. Only support mode 2 for now.
  if (buffer.i[3] != 2) {
    mprinterr("Error: Mode %i; currently only mode 2 for CCP4 files is supported.\n", buffer.i[3]);
    return 1;
  }
  // Check offsets.
  if (buffer.i[4] != 0 || buffer.i[5] != 0 || buffer.i[6] != 0)
    mprintf("Warning: Non-zero offsets present. This is not yet supported and will be ignored.\n");
  // Check that mapping is col=x row=y section=z
  if (buffer.i[16] != 1 || buffer.i[17] != 2 || buffer.i[18] != 3) {
    mprinterr("Error: Currently only support cols=X, rows=Y, sections=Z\n");
    return 1;
  }
  if (buffer.i[24] != 0) {
    mprintf("Warning: Skew information present but not yet supported and will be ignored.\n");
    return 1;
  }

  // Read 10 80 character text labels
  char Labels[801];
  Labels[800] = '\0';
  infile.Read( Labels, 200*wSize );
  mprintf("\t%s\n", Labels);
  // Symmetry records: operators separated by * and grouped into 'lines' of 80 characters
  int NsymmRecords = buffer.i[23] / 80;
  if (NsymmRecords > 0) {
    char symBuffer[81];
    mprintf("\t%i symmetry records.\n", NsymmRecords);
    for (int ib = 0; ib != NsymmRecords; ib++) {
      infile.Gets( symBuffer, 80 );
      mprintf("\t%s\n", symBuffer);
    }
  }

  // Add grid data set. Default to float for now.
  DataSet* gridDS = datasetlist.AddSet( DataSet::GRID_FLT, dsname, "GRID" );
  if (gridDS == 0) return 1;
  DataSet_GridFlt& grid = static_cast<DataSet_GridFlt&>( *gridDS );
  // Allocate grid from dims and spacing. FIXME OK to assume zero origin?
  if (grid.Allocate_N_O_Box( buffer.i[7], buffer.i[8], buffer.i[9],
                             Vec3(0.0), Box(buffer.f + 10) ) != 0)
  {
    mprinterr("Error: Could not allocate grid.\n");
    return 1;
  }
  // FIXME: Grids are currently indexed so Z is fastest changing.
  //        Should be able to change indexing in grid DataSet.
  size_t mapSize = buffer.i[7] * buffer.i[8] * buffer.i[9];
  mprintf("\tCCP4 map has %zu elements\n", mapSize);
  mprintf("\tDensity: Min=%f  Max=%f  Mean=%f  RMS=%f\n",
          buffer.f[19], buffer.f[20], buffer.f[21], buffer.f[54]);
  std::vector<float> mapbuffer( mapSize );
  int mapBytes = mapSize * wSize;
  int numRead = infile.Read( &mapbuffer[0], mapBytes );
  if (numRead < 1) {
    mprinterr("Error: Could not read CCP4 map data.\n");
    return 1;
  } else if (numRead < mapBytes)
    mprintf("Warning: Expected %i bytes, read only %i bytes\n", mapBytes, numRead);
  if (isBigEndian) endian_swap(&mapbuffer[0], mapSize);

  // FIXME: Place data into grid DataSet with correct ordering.
  int gidx = 0;
  int NXY = buffer.i[7] * buffer.i[8];
  for (int ix = 0; ix != buffer.i[7]; ix++)
    for (int iy = 0; iy != buffer.i[8]; iy++)
      for (int iz = 0; iz != buffer.i[9]; iz++) {
        int midx = (iz * NXY) + (iy * buffer.i[7]) + ix;
        grid[gidx++] = mapbuffer[midx];
      }

  infile.CloseFile();
  return 0;
}

// DataIO_CCP4::WriteHelp()
void DataIO_CCP4::WriteHelp() {
  mprintf("\ttitle <title>: Set CCP4 output title.\n");
}

// DataIO_CCP4::processWriteArgs()
int DataIO_CCP4::processWriteArgs(ArgList& argIn) {
  title_ = argIn.GetStringKey("title");
  return 0;
}

// DataIO_CCP4::WriteData()
int DataIO_CCP4::WriteData(FileName const& fname, DataSetList const& setList)
{
  // Open output file
  CpptrajFile outfile;
  if (outfile.OpenWrite(fname)) {
    mprinterr("Error: Could not open CCP4 output file '%s'.\n", fname.full());
    return 1;
  }
  // Warn about writing multiple sets
  if (setList.size() > 1)
    mprintf("Warning: %s: Writing multiple 3D sets in CCP4 format not supported.\n"
            "Warning:   Only writing first set.\n", fname.full());
  return WriteSet3D( setList.begin(), outfile );
}

// DataIO_CCP4::WriteSet3D()
int DataIO_CCP4::WriteSet3D( DataSetList::const_iterator const& setIn, CpptrajFile& outfile ) {
  if ((*setIn)->Size() < 1) return 1; // SANITY CHECK: No empty grid allowed
  if ((*setIn)->Ndim() != 3) {
    mprinterr("Internal Error: DataSet %s in DataFile %s has %zu dimensions, expected 3.\n",
              (*setIn)->legend(), outfile.Filename().full(), (*setIn)->Ndim());
    return 1;
  }
  DataSet_3D const& grid = static_cast<DataSet_3D const&>( *(*setIn) );
  // Check input grid
  Vec3 OXYZ = grid.Bin().GridOrigin();
  if (OXYZ[0] < 0.0 || OXYZ[1] < 0.0 || OXYZ[2] < 0.0 ||
      OXYZ[0] > 0.0 || OXYZ[1] > 0.0 || OXYZ[2] > 0.0)
    mprintf("Warning: Grid '%s' origin is not 0.0, 0.0, 0.0\n"
            "Warning:  Origin other than 0.0 not yet supported for CCP4 write.\n");
  // Set default title if none set 
  if (title_.empty())
    title_.assign("CPPTRAJ CCP4 map volumetric data, set '" + grid.Meta().Legend() +
                  "'. Format revision A.");
  // Check that title is not too big
  if (title_.size() > 800) {
    mprintf("Warning: CCP4 title is too large, truncating.\n");
    title_.resize( 800 );
  }
  // Set up and write header
  headerbyte buffer;
  buffer.i[0] = (int)grid.NX();
  buffer.i[1] = (int)grid.NY();
  buffer.i[2] = (int)grid.NZ();
  buffer.i[3] = 2; // Only mode 2 supported
  buffer.i[4] = 0; // No offsets
  buffer.i[5] = 0;
  buffer.i[6] = 0;
  buffer.i[7] = (int)grid.NX();
  buffer.i[8] = (int)grid.NY();
  buffer.i[9] = (int)grid.NZ();
  Box box( grid.Bin().Ucell() );
  buffer.f[10] = (float)box[0];
  buffer.f[11] = (float)box[1];
  buffer.f[12] = (float)box[2];
  buffer.f[13] = (float)box[3];
  buffer.f[14] = (float)box[4];
  buffer.f[15] = (float)box[5];
  buffer.i[16] = 1; // Cols = X
  buffer.i[17] = 2; // Rows = Y
  buffer.i[18] = 3; // Secs = Z
  // Determine min, max, and mean of data
  double mean = grid[0];
  double gmin = grid[0];
  double gmax = grid[0];
  double rmsd = grid[0] * grid[0];
  for (unsigned int i = 1; i < grid.Size(); i++) {
    gmin = std::min(grid[i], gmin);
    gmax = std::max(grid[i], gmax);
    mean += grid[i];
    rmsd += grid[i] * grid[i];
  }
  mean /= (double)grid.Size();
  rmsd /= (double)grid.Size();
  rmsd = rmsd - (mean * mean);
  if (rmsd > 0.0)
    rmsd = sqrt(rmsd);
  else
    rmsd = 0.0;
  mprintf("\t%s\n", title_.c_str());
  mprintf("\tDensity: Min=%f  Max=%f  Mean=%f  RMS=%f\n", gmin, gmax, mean, rmsd);
  buffer.f[19] = (float)gmin;
  buffer.f[20] = (float)gmax;
  buffer.f[21] = (float)mean;
  buffer.i[22] = 1; // Assume P1
  buffer.i[23] = 0; // No bytes for symmetry ops
  buffer.i[24] = 0; // No skew transform
  // Skew matrix (S11, S12, S13, S21, ...) and translation; 12 total, followed
  // by 15 'future use'; zero all.
  std::fill( buffer.i+25, buffer.i+52, 0 );
  // MAP and machine precision. FIXME determine endianness!
  buffer.c[208] = 'M';
  buffer.c[209] = 'A';
  buffer.c[210] = 'P';
  buffer.c[211] = ' ';
  buffer.c[212] = 0x44; // little endian
  buffer.c[213] = 0x41;
  buffer.c[214] = 0x00;
  buffer.c[215] = 0x00;
  // Determine RMS deviation from mean.
  buffer.f[54] = (float)rmsd;
  // Determine number of labels being used
  buffer.i[55] = (int)title_.size() / 80;
  if ( ((int)title_.size() % 80) != 0)
    buffer.i[55]++;
  outfile.Write( buffer.c, 224*sizeof(unsigned char) );

  // Write labels; 10 lines, 80 chars each
  outfile.Write( title_.c_str(), title_.size() );
  // FIXME this seems wasteful.
  std::vector<char> remainder( 800 - title_.size(), 0 );
  outfile.Write( &remainder[0], remainder.size() );
  remainder.clear();

  // No symmetry bytes

  // Store data in buffer, then write. X changes fastest.
  // SANITY CHECK; This will result in invalid files if size of float is not 4.
  if (sizeof(float) != wSize)
    mprintf("Warning: Size of float on this system is %zu, not 4.\n"
            "Warning:  Resulting CCP4 file data will not conform to standard.\n", sizeof(float));
  std::vector<float> mapbuffer( grid.Size() );
  std::vector<float>::iterator it = mapbuffer.begin();
  for (unsigned int iz = 0; iz != grid.NZ(); iz++)
    for (unsigned int iy = 0; iy != grid.NY(); iy++)
      for (unsigned int ix = 0; ix != grid.NX(); ix++)
        *(it++) = grid.GetElement( ix, iy, iz );
  outfile.Write( &mapbuffer[0], mapbuffer.size() * sizeof(float) );

  outfile.CloseFile();
  return 0;
}
