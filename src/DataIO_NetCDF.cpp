#include "DataIO_NetCDF.h"
#include "CpptrajStdio.h"
#ifdef BINTRAJ
# include <map>
# include <string>
# include <netcdf.h>
# include "NC_Routines.h"
# include "Version.h"
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
/*class NC_dimension {
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
};*/

/// Hold a pool of pointers to DataSets in the list.
class DataIO_NetCDF::SetPool {
  public:
    /// CONSTRUCTOR - place sets from DataSetList in this pool
    SetPool(DataSetList const& dsl) : nUsed_(0) {
      sets_.reserve( dsl.size() );
      isUsed_.assign( dsl.size(), false );
      for (DataSetList::const_iterator it = dsl.begin(); it != dsl.end(); ++it)
        sets_.push_back( *it );
    }
    /// \return set at idx
    DataSet const* Set(unsigned int idx) const { return sets_[idx]; }
    /// \return Number of sets
    unsigned int Nsets() const { return sets_.size(); }
    /// \return true if set at idx has been used
    bool IsUsed(unsigned int idx) const { return isUsed_[idx]; }
    /// \return true if all sets have been marked as used
    bool AllUsed() const { return (nUsed_ == isUsed_.size()); }

    /// Mark set at idx as used
    void MarkUsed(unsigned int idx) {
      isUsed_[idx] = true;
      nUsed_++;
    }
  private:
    std::vector<DataSet const*> sets_;
    std::vector<bool> isUsed_;
    unsigned int nUsed_;
};

// DataIO_NetCDF::WriteData()
int DataIO_NetCDF::WriteData(FileName const& fname, DataSetList const& dsl)
{
# ifdef BINTRAJ
  int ncid = -1;
  // TODO check existing file
  if (NC::CheckErr( nc_create( fname.full(), NC_64BIT_OFFSET, &ncid ) ))
    return 1;

  // Place incoming DataSets into a pool. As they are handled they will
  // be removed from the pool.
  SetPool setPool( dsl );

  typedef std::vector<DataSet const*> SetArray;

  // Check our incoming data sets. Try to find common dimensions.
  for (unsigned int idx = 0; idx < setPool.Nsets(); idx++)
  {
    if (setPool.IsUsed(idx)) continue;

    DataSet const* ds = setPool.Set( idx );

    if (ds->Group() == DataSet::SCALAR_1D) {
      // ----- 1D scalar -------------------------
      SetArray sets(1, ds);
      setPool.MarkUsed( idx );
      // Group this set with others that share the same dimension and length.
      Dimension const& dim = ds->Dim(0);
      for (unsigned int jdx = idx + 1; jdx < setPool.Nsets(); jdx++)
      {
        DataSet const* ds2 = setPool.Set( jdx );
        if (ds->Size() == ds2->Size() && dim == ds2->Dim(0)) {
          sets.push_back( ds2 );
          setPool.MarkUsed( jdx );
        }
      }
      mprintf("DEBUG: Sets for dimension '%s' %f %f:", dim.label(), dim.Min(), dim.Step());
      for (SetArray::const_iterator it = sets.begin(); it != sets.end(); ++it)
        mprintf(" %s", (*it)->legend());
      mprintf("\n");
    } else {
      mprinterr("Error: '%s' is an unhandled set type for NetCDF.\n", ds->legend());
      return 1;
    }
  } // END loop over set pool
  if (!setPool.AllUsed()) {
    mprintf("Warning: Not all sets were used.\n");
  }
/*
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
*/
  // Attributes
  //if (NC::CheckErr(nc_put_att_text(ncid_,NC_GLOBAL,"title",title.size(),title.c_str())) ) {
  //  mprinterr("Error: Writing title.\n");
  //  return 1;
  //}
  if (NC::CheckErr(nc_put_att_text(ncid,NC_GLOBAL,"application",5,"AMBER")) ) {
    mprinterr("Error: Writing application.\n");
    return 1;
  }
  if (NC::CheckErr(nc_put_att_text(ncid,NC_GLOBAL,"program",7,"cpptraj")) ) {
    mprinterr("Error: Writing program.\n");
    return 1;
  }
  std::string programVersion(CPPTRAJ_INTERNAL_VERSION);
  if (NC::CheckErr(nc_put_att_text(ncid,NC_GLOBAL,"programVersion",
                                   programVersion.size(), programVersion.c_str())))
  {
    mprinterr("Error: Writing program version.\n");
    return 1;
  }

  // End netcdf definitions
  if (NC::CheckErr(nc_enddef(ncid))) {
    mprinterr("NetCDF data error on ending definitions.");
    return 1;
  }
  NC::Debug(ncid);
  return 0;
# else
  return 1;
# endif
}
