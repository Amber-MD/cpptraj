#include "Exec_SortEnsembleData.h"
#include "CpptrajStdio.h"

// Exec_SortEnsembleData::Help()
void Exec_SortEnsembleData::Help() const
{

}

inline bool CheckError(int err) {
# ifdef MPI
  if (Parallel::EnsembleComm().CheckError( err )) return 1;
# else
  if (err != 0) return 1;
# endif
  return 0;
}

//  Exec_SortEnsembleData::Sort_pH_Data()
int Exec_SortEnsembleData::Sort_pH_Data(DataSetList const& setsToSort) const {
  return 0;
}

// Exec_SortEnsembleData::SortData()
int Exec_SortEnsembleData::SortData(DataSetList const& setsToSort) const {
  int err = 0;
  if (setsToSort.empty()) {
    rprinterr("Error: No sets selected.\n");
    err = 1;
  }
  if (CheckError(err)) return 1;

  DataSet::DataType dtype = setsToSort[0]->Type();
  for (DataSetList::const_iterator ds = setsToSort.begin(); ds != setsToSort.end(); ++ds)
    if (dtype != (*ds)->Type()) {
      rprinterr("Error: Set '%s' has different type than first set.\n", (*ds)->legend());
      err = 1;
      break;
    }
  if (CheckError(err)) return 1; 

# ifdef MPI
  typedef std::vector<int> Iarray;
  Iarray Dtypes( Parallel::EnsembleComm().Size(), -1 );
  if ( Parallel::EnsembleComm().AllGather( &dtype, 1, MPI_INT, &Dtypes[0] ) ) return 1;
  for (int rank = 1; rank < Parallel::EnsembleComm().Size(); rank++)
    if (Dtypes[0] != Dtypes[rank]) {
      mprinterr("Error: Set types on rank %i do not match types on rank 0.\n", rank);
      err = 1;
      break;
    }
  if (Parallel::EnsembleComm().CheckError( err )) return 1;
# endif

  // Only work for pH data for now.
  if (dtype != DataSet::PH) {
    mprinterr("Error: Only works for pH data for now.\n");
    return 1;
  }

  err = Sort_pH_Data( setsToSort );

  return err;
}

// Exec_SortEnsembleData::Execute()
Exec::RetType Exec_SortEnsembleData::Execute(CpptrajState& State, ArgList& argIn)
{
  rprintf("DEBUG: Entering sortensembledata.\n");
  DataSetList setsToSort;
  std::string dsarg = argIn.GetStringNext();
  while (!dsarg.empty()) {
    setsToSort += State.DSL().GetMultipleSets( dsarg );
    dsarg = argIn.GetStringNext();
  }
  setsToSort.List();

  int err = 0;
# ifdef MPI
  // For now, require ensemble mode in parallel.
  if (Parallel::EnsembleComm().IsNull()) {
    mprinterr("Error: Data set ensemble sort requires ensemble mode in parallel.\n");
    return CpptrajState::ERR;
  }
  // If not a TrajComm master we do not have complete data, so exit now.
  if (Parallel::TrajComm().Master())
    err = SortData( setsToSort );
  if (Parallel::World().CheckError( err ))
# else
  if (SortData( setsToSort )) 
# endif
    return CpptrajState::ERR;
  return CpptrajState::OK;
}
