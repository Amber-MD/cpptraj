#include "Exec_SortEnsembleData.h"
#include "CpptrajStdio.h"
#include "DataSet_PH.h"

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
  // Gather pH data
  typedef std::vector<DataSet_PH*> Parray;
  Parray PH;
  for (DataSetList::const_iterator ds = setsToSort.begin(); ds != setsToSort.end(); ++ds)
    PH.push_back( (DataSet_PH*)*ds );

  typedef std::vector<double> Darray;
  Darray pHvalues;
# ifdef MPI
  pHvalues.resize( Parallel::Ensemble_Size() );
  Darray phtmp;
  for (Parray::const_iterator ds = PH.begin(); ds != PH.end(); ++ds)
    phtmp.push_back( (*ds)->pH_Values()[0] );
  if (Parallel::EnsembleComm().AllGather(&phtmp[0], phtmp.size(), MPI_DOUBLE, &pHvalues[0])) {
    rprinterr("Error: Gathering pH values.\n");
    return 1;
  }
# else
  for (Parray::const_iterator ds = PH.begin(); ds != PH.end(); ++ds)
    pHvalues.push_back( (*ds)->pH_Values()[0] );
# endif
  ReplicaInfo::Map<double> pH_map;
  if (pH_map.CreateMap( pHvalues )) {
    rprinterr("Error: Duplicate pH value detected (%.2f) in ensemble.\n", pH_map.Duplicate());
    return 1;
  }
  mprintf("\tInitial pH values:");
  for (ReplicaInfo::Map<double>::const_iterator ph = pH_map.begin(); ph != pH_map.end(); ++ph)
    mprintf(" %6.2f", *ph);
  mprintf("\n");
  
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
# ifdef MPI
  // Number of sets to sort should be equal to # members I am responsible for.
  if (Parallel::N_Ens_Members() != (int)setsToSort.size()) {
    rprinterr("Internal Error: Number of ensemble members (%i) != # sets to sort (%zu)\n",
               Parallel::N_Ens_Members(), setsToSort.size());
    return 1;
  }
# endif

  DataSet::DataType dtype = setsToSort[0]->Type();
  for (DataSetList::const_iterator ds = setsToSort.begin(); ds != setsToSort.end(); ++ds) {
    rprintf("\t%s\n", (*ds)->legend());
    if ((*ds)->Size() < 1) {
      rprinterr("Error: Set '%s' is empty.\n", (*ds)->legend());
      err = 1;
      break;
    }
    if (dtype != (*ds)->Type()) {
      rprinterr("Error: Set '%s' has different type than first set.\n", (*ds)->legend());
      err = 1;
      break;
    }
  }
  if (CheckError(err)) return 1; 

# ifdef MPI
  Parallel::EnsembleComm().Barrier(); // DEBUG
  typedef std::vector<int> Iarray;
  Iarray Dtypes( Parallel::EnsembleComm().Size(), -1 );
  if ( Parallel::EnsembleComm().AllGather( &dtype, 1, MPI_INT, &Dtypes[0] ) ) return 1;
  for (int rank = 1; rank < Parallel::EnsembleComm().Size(); rank++)
    if (Dtypes[0] != Dtypes[rank]) {
      rprinterr("Error: Set types on rank %i do not match types on rank 0.\n", rank);
      err = 1;
      break;
    }
  if (Parallel::EnsembleComm().CheckError( err )) return 1;
# endif

  // Only work for pH data for now.
  if (dtype != DataSet::PH) {
    rprinterr("Error: Only works for pH data for now.\n");
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
    rprinterr("Error: Data set ensemble sort requires ensemble mode in parallel.\n");
    return CpptrajState::ERR;
  }
  // Only TrajComm masters have complete data.
  if (Parallel::TrajComm().Master())
    err = SortData( setsToSort );
  if (Parallel::World().CheckError( err ))
# else
  if (SortData( setsToSort )) 
# endif
    return CpptrajState::ERR;
  return CpptrajState::OK;
}
