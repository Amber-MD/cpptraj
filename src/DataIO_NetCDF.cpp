#include "DataIO_NetCDF.h"
#include "CpptrajStdio.h"
#ifdef BINTRAJ
# include <set>
# include <netcdf.h>
# include "NC_Routines.h"
#endif

/// CONSTRUCTOR
DataIO_NetCDF::DataIO_NetCDF() :
  DataIO(true, true, true) // Valid for 1D, 2D, 3D
{

}

// DataIO_NetCDF::ID_DataFormat()
bool DataIO_NetCDF::ID_DataFormat(CpptrajFile& infile)
{
  if (infile.OpenFile()) return false;
  // Read first 3 bytes
  unsigned char magic[3];
  magic[0] = 0;
  magic[1] = 0;
  magic[2] = 0;
  infile.Read(magic, 3);
  infile.CloseFile();
  if (magic[0] == 'C' && magic[1] == 'D' && magic[2] == 'F') {
#   ifdef BINTRAJ
    return true;
#   else
    mprintf("Warning: '%s' is a NetCDF file but CPPTRAJ was compiled without NetCDF support.\n",
            infile.Filename().full());
#   endif
  }
  return false;
}

// DataIO_NetCDF::ReadHelp()
void DataIO_NetCDF::ReadHelp()
{

}

// DataIO_NetCDF::processReadArgs()
int DataIO_NetCDF::processReadArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_NetCDF::ReadData()
int DataIO_NetCDF::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{

  return 1;
}

// DataIO_NetCDF::WriteHelp()
void DataIO_NetCDF::WriteHelp()
{

}

// DataIO_NetCDF::processWriteArgs()
int DataIO_NetCDF::processWriteArgs(ArgList& argIn)
{

  return 0;
}

class NC_dimension {
  public:
    NC_dimension(std::string const& l, int s) : label_(l), size_(s) {}

    bool operator==(NC_dimension const& rhs) const {
      if (rhs.size_ != size_) return false;
      if (rhs.label_ != label_) return false;
      return true;
    }

    bool operator<(NC_dimension const& rhs) const {
      if (rhs.size_ == size_) {
        return (label_ < rhs.label_);
      } else {
        return (size_ < rhs.size_);
      }
    }

  std::string const& Label() const { return label_; }
  int Size() const { return size_; }

  private:

    std::string label_;
    int size_;
};

// DataIO_NetCDF::WriteData()
int DataIO_NetCDF::WriteData(FileName const& fname, DataSetList const& dsl)
{
# ifdef BINTRAJ
  int ncid = -1;
  // TODO check existing file
  if (NC::CheckErr( nc_create( fname.full(), NC_64BIT_OFFSET, &ncid ) ))
    return 1;

  // Check our incoming data sets. Try to find common dimensions.
  std::set<NC_dimension> dims_;

  for (DataSetList::const_iterator dsit = dsl.begin(); dsit != dsl.end(); ++dsit)
  {
    DataSet const* ds = *dsit;
    if (ds->Group() == DataSet::SCALAR_1D) {
      // 1 dimension
      Dimension const& dim = ds->Dim(0);

      dims_.insert( NC_dimension(dim.Label(), ds->Size()) );
    }
  }
  mprintf("DEBUG: Dimensions:\n");
  for (std::set<NC_dimension>::const_iterator it = dims_.begin(); it != dims_.end(); ++it)
    mprintf("\t'%s' %i\n", it->Label().c_str(), it->Size());

  return 0;
# else
  return 1;
# endif
}
