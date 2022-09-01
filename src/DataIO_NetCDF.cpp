#include "DataIO_NetCDF.h"
#include "CpptrajStdio.h"
#ifdef BINTRAJ
# include <map>
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

/// Used to track unique DataSet dimensions
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
  std::vector<DataSet const*> const& Sets() const { return sets_; }

  void AddSet(DataSet const* ds) { sets_.push_back( ds ); }
  private:

    std::string label_;
    int size_;
    std::vector<DataSet const*> sets_;
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
  // Give each unique dimension an index.
  typedef std::pair<NC_dimension,int> DimPair;
  typedef std::map<NC_dimension,int> DimSet;
  DimSet dims_;
  // Hold dims for each set
  typedef std::vector<int> Iarray;
  typedef std::vector<Iarray> DimArray;
  DimArray dimIndices_;

  for (DataSetList::const_iterator dsit = dsl.begin(); dsit != dsl.end(); ++dsit)
  {
    DataSet const* ds = *dsit;
    if (ds->Group() == DataSet::SCALAR_1D) {
      // 1 dimension
      Dimension const& dim = ds->Dim(0);
      NC_dimension ncdim(dim.Label(), ds->Size());

      DimSet::iterator ret = dims_.find( ncdim );

      if (ret == dims_.end()) {
        std::pair<DimSet::iterator,bool> it = dims_.insert( DimPair(ncdim,dims_.size()) );
        ret = it.first;
      }
      dimIndices_.push_back( Iarray(1, ret->second) );

      //std::pair<DimSet::iterator,bool> ret = dims_.insert( NC_dimension(dim.Label(), ds->Size()) );
      //mprintf("%i\n", ret.first->Size());
      //ret.first->AddSet( ds );
    }
  }
  mprintf("DEBUG: Dimensions:\n");
  for (DimSet::const_iterator it = dims_.begin(); it != dims_.end(); ++it)
    mprintf("\t%i : '%s' %i\n", it->second, it->first.Label().c_str(), it->first.Size());
  mprintf("DEBUG: Set dim indices:\n");
  DimArray::const_iterator didx = dimIndices_.begin();
  for (DataSetList::const_iterator dsit = dsl.begin(); dsit != dsl.end(); ++dsit, ++didx)
  {
    mprintf("\t%s :", (*dsit)->legend());
    for (Iarray::const_iterator it = didx->begin(); it != didx->end(); ++it)
      mprintf(" %i", *it);
    mprintf("\n");
  }

  // Define dimensions
  Iarray dimIds;
  dimIds.resize( dims_.size() );
  int* dimIdsPtr = &dimIds[0];
  for (DimSet::const_iterator it = dims_.begin(); it != dims_.end(); ++it, ++dimIdsPtr) {
    if (NC::CheckErr( nc_def_dim(ncid, it->first.Label().c_str(), it->first.Size(), dimIdsPtr ))) {
      mprinterr("Error: Could not define dimension '%s'\n", it->first.Label().c_str());
      return 1;
    }
  }

  return 0;
# else
  return 1;
# endif
}
