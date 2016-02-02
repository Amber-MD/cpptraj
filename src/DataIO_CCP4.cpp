#include "DataIO_CCP4.h"
#include "CpptrajStdio.h"
#include "ByteRoutines.h"

const size_t DataIO_CCP4::wSize = 4;

bool DataIO_CCP4::MapCharsValid(const unsigned char* MAP) {
  return (MAP[0] == 'M' && MAP[1] == 'A' && MAP[2] == 'P' && MAP[3] == ' ');
}

// DataIO_CCP4::ID_DataFormat()
bool DataIO_CCP4::ID_DataFormat( CpptrajFile& infile ) {
  bool isCCP4 = false;
  if (!infile.OpenFile()) {
    unsigned char MAP[4];
    if (infile.Seek(52 * wSize) == 0) {
      infile.Read( MAP, wSize );
      isCCP4 = MapCharsValid( MAP );
    }
    infile.CloseFile();
  }
  return isCCP4;
}

//TODO may need to use byte swapping routines
// DataIO_CCP4::ReadData()
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
    mprintf("DEBUG: Little endian.\n");
    // SANITY CHECK
    if ( !(buffer.c[212] == 0x44 && buffer.c[213] == 0x41 &&
           buffer.c[214] == 0x00 && buffer.c[215] == 0x00) )
      mprintf("Warning: Invalid machine stamp: %x %x %x %x : assuming little endian.\n",
              buffer.c[212], buffer.c[213], buffer.c[214], buffer.c[215]);
  } else {
    mprintf("DEBUG: Big endian.\n");
    // Perform endian swapping on header if necessary
    endian_swap(buffer.i, 56);
  }

  // Read columns, rows, sections (fastest changing to slowest)
  mprintf("DEBUG: Columns=%i  Rows=%i  Sections=%i\n", buffer.i[0], buffer.i[1], buffer.i[2]);
  // Read mode: 0 = envelope stored as signed bytes (from -128 lowest to 127 highest)
  //            1 = Image     stored as Integer*2
  //            2 = Image     stored as Reals
  //            3 = Transform stored as Complex Integer*2
  //            4 = Transform stored as Complex Reals
  //            5 == 0
  // Only support mode 2 for now.
  mprintf("DEBUG: Mode=%i\n", buffer.i[3]);
  if (buffer.i[3] != 2) {
    mprinterr("Error: Mode %i; currently only mode 2 for CCP4 files is supported.\n", buffer.i[3]);
    return 1;
  }
  // Column, row, and section offsets
  mprintf("DEBUG: Offsets: C=%i  R=%i  S=%i\n", buffer.i[4], buffer.i[5], buffer.i[6]);
  // TODO check non-zero offsets?
  // Intervals along X, Y, Z
  mprintf("DEBUG: NXYZ={ %i %i %i }\n", buffer.i[7], buffer.i[8], buffer.i[9]);
  // Read 3x cell lengths (Ang) and 3x cell angles (deg)
  mprintf("DEBUG: Box XYZ={ %f %f %f }  ABG={ %f %f %f }\n",
          buffer.f[10], buffer.f[11], buffer.f[12],
          buffer.f[13], buffer.f[14], buffer.f[15]);
  // Read map of which axes correspond to cols, rows, sections (1,2,3 = x,y,z)
  mprintf("DEBUG: Map: ColAxis=%i  RowAxis=%i  SecAxis=%i\n",
          buffer.i[16], buffer.i[17], buffer.i[18]);
  // TODO only allow 1,2,3
  // Read min, max, and mean density
  mprintf("DEBUG: Density: Min=%f  Max=%f  Mean=%f\n", buffer.f[19], buffer.f[20], buffer.f[21]);
  // Read space group, bytes used for storing symm ops, flag for skew transform
  mprintf("DEBUG: SpaceGroup#=%i  SymmOpBytes=%i  SkewFlag=%i\n",
          buffer.i[22], buffer.i[23], buffer.i[24]);
  // If skew flag != 0, skew transformation is from standard orthogonal
  // coordinate frame (as used for atoms) to orthogonal map frame, as:
  //     Xo(map) = S * (Xo(atoms) - t)
  // Read skew matrix (in order S11, S12, S13, S21 etc)
  const int* MSKEW = buffer.i + 25;
  mprintf("DEBUG: Skew matrix: %i %i %i\n"
          "                    %i %i %i\n"
          "                    %i %i %i\n", MSKEW[0], MSKEW[1], MSKEW[2], MSKEW[3],
          MSKEW[4], MSKEW[5], MSKEW[6], MSKEW[7], MSKEW[8]);
  // Read skew translation.
  const int* TSKEW = buffer.i + 34;
  mprintf("DEBUG: Skew translation: %i %i %i\n", TSKEW[0], TSKEW[1], TSKEW[2]);
  // The next 15 values are for future use and can be skipped.
  // Read RMS deviation of map from mean and # labels
  mprintf("DEBUG: RMSD=%f  Nlabels=%i\n", buffer.f[54], buffer.i[55]);
  // 10 80 character text labels
  char Labels[801];
  Labels[800] = '\0';
  infile.Read( Labels, 200*wSize );
  mprintf("DEBUG: Labels:\n%s\n", Labels);
  // Symmetry records: operators separated by * and grouped into 'lines' of 80 characters
  int NsymmRecords = buffer.i[23] / 80;
  char symBuffer[81];
  mprintf("DEBUG: %i symmetry records.\n", NsymmRecords);
  for (int ib = 0; ib != NsymmRecords; ib++) {
    infile.Gets( symBuffer, 80 );
    mprintf("\t%s\n", symBuffer);
  }

  infile.CloseFile();
  return 0;
}

// DataIO_CCP4::WriteHelp()
void DataIO_CCP4::WriteHelp() {

}

// DataIO_CCP4::processWriteArgs()
int DataIO_CCP4::processWriteArgs(ArgList& argIn) {
  return 0;
}

// DataIO_CCP4::WriteData()
int DataIO_CCP4::WriteData(FileName const& fname, DataSetList const& setList)
{
  return 1;
}
