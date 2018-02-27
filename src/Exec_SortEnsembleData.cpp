#include "Exec_SortEnsembleData.h"
#include "CpptrajStdio.h"
#include "DataSet_pH_REMD.h"
#include "DataSet_pH.h"
#include "StringRoutines.h" // doubleToString

inline bool CheckError(int err) {
# ifdef MPI
  if (Parallel::EnsembleComm().CheckError( err )) return 1;
# else
  if (err != 0) return 1;
# endif
  return 0;
}

//  Exec_SortEnsembleData::Sort_pH_Data()
int Exec_SortEnsembleData::Sort_pH_Data(DataSetList const& setsToSort, DataSetList& OutputSets,
                                        unsigned int maxFrames)
const
{
  // Cast sets back to DataSet_pH_REMD
  typedef std::vector<DataSet_pH_REMD*> Parray;
  Parray PHsets;
  for (DataSetList::const_iterator ds = setsToSort.begin(); ds != setsToSort.end(); ++ds)
    PHsets.push_back( (DataSet_pH_REMD*)*ds );

  // Gather initial pH data values, ensure no duplicates
  typedef std::vector<double> Darray;
  Darray pHvalues;
# ifdef MPI
  pHvalues.resize( Parallel::Ensemble_Size() );
  Darray phtmp;
  for (Parray::const_iterator ds = PHsets.begin(); ds != PHsets.end(); ++ds)
    phtmp.push_back( (*ds)->pH_Values()[0] );
  if (Parallel::EnsembleComm().AllGather(&phtmp[0], phtmp.size(), MPI_DOUBLE, &pHvalues[0])) {
    rprinterr("Error: Gathering pH values.\n");
    return 1;
  }
# else
  for (Parray::const_iterator ds = PHsets.begin(); ds != PHsets.end(); ++ds)
    pHvalues.push_back( (*ds)->pH_Values()[0] );
# endif
  ReplicaInfo::Map<double> pH_map;
  if (pH_map.CreateMap( pHvalues )) {
    rprinterr("Error: Duplicate pH value detected (%.2f) in ensemble.\n", pH_map.Duplicate());
    return 1;
  }
  Darray sortedPH;
  mprintf("\tInitial pH values:");
  for (ReplicaInfo::Map<double>::const_iterator ph = pH_map.begin(); ph != pH_map.end(); ++ph)
  {
    mprintf(" %6.2f", ph->first);
    sortedPH.push_back( ph->first );
  }
  mprintf("\n");

  // Create sets to hold sorted pH values. Create a set for each pH value
  // and each residue. Final output sets will be PH0R0, PH0R1, PH1R0, ...
  // TODO check that residue info all the same
  DataSet_pH_REMD::Rarray const& Residues = PHsets[0]->Residues();
  if (debug_ > 0)
    rprintf("DEBUG: Sorting %u frames for %zu sets, %zu pH values.\n",
            maxFrames, PHsets.size(), pHvalues.size());
  for (unsigned int idx = 0; idx != sortedPH.size(); idx++) {
    OutputSets.SetEnsembleNum( idx );
    for (unsigned int res = 0; res != Residues.size(); ++res) {
      MetaData md(PHsets[0]->Meta().Name(), Residues[res].Name().Truncated(), Residues[res].Num());
      DataSet_pH* out = (DataSet_pH*)OutputSets.AddSet( DataSet::PH, md );
      if (out==0) return 1;
      //out->SetLegend( "pH " + doubleToString( sortedPH[idx] ) );
      out->Set_Solvent_pH( sortedPH[idx] );
      out->SetResidueInfo( Residues[res] );
      out->SetTimeValues(PHsets[0]->Time());
      out->Resize( maxFrames );
    }
  }

  // Loop over unsorted sets
  for (Parray::const_iterator ds = PHsets.begin(); ds != PHsets.end(); ++ds)
  {
    unsigned int phidx = 0;
    for (unsigned int n = 0; n < maxFrames; n++)
    {
      float phval = (*ds)->pH_Values()[n];
      int setidx = pH_map.FindIndex( phval ) * Residues.size();
      //rprintf("DEBUG: %6u Set %10s pH= %6.2f going to %2i\n", n+1, (*ds)->legend(), phval, idx);
      //mflush();
      for (unsigned int res = 0; res < (*ds)->Residues().size(); res++, setidx++, phidx++)
      {
        DataSet_pH* out = (DataSet_pH*)OutputSets[setidx];
        //if (res == 0 && idx == 0) {
        //  rprintf("DEBUG: Frame %3u res %2u State %2i pH %6.2f\n", 
        //          n, res, (*ds)->Res(res).State(n), phval);
        //  mflush();
        //}
        out->SetState(n, (*ds)->ResStates()[phidx], (*ds)->RecordType(n));
      }
    }
  }
# ifdef MPI
  // Now we need to reduce down each set onto the thread where it belongs.
  if (Parallel::World().Size() > 1) {
    for (int idx = 0; idx != (int)OutputSets.size(); idx++) {
      DataSet_pH* out = (DataSet_pH*)OutputSets[idx];
      int ensembleRank = Parallel::MemberEnsCommRank( out->Meta().EnsembleNum() );
      //rprintf("DEBUG: Consolidate set %s to rank %i\n", out->legend(), ensembleRank);
      out->Consolidate( Parallel::EnsembleComm(), ensembleRank );
    }
    // Remove sets that do not belong on this rank
    for (int idx = (int)OutputSets.size() - 1; idx > -1; idx--) {
      DataSet* out = OutputSets[idx];
      int ensembleRank = Parallel::MemberEnsCommRank( out->Meta().EnsembleNum() );
      if (ensembleRank != Parallel::EnsembleComm().Rank()) {
        //rprintf("DEBUG: Remove set %s (%i) from rank %i\n", out->legend(),
        //        idx, Parallel::EnsembleComm().Rank());
        OutputSets.RemoveSet( out );
      }
    }
  }
# endif
  return 0;
}

// Exec_SortEnsembleData::SortData()
int Exec_SortEnsembleData::SortData(DataSetList const& setsToSort, DataSetList& OutputSets)
const
{
  int err = 0;
  if (setsToSort.empty()) {
    rprinterr("Error: No sets selected.\n");
    err = 1;
  }
  if (CheckError(err)) return 1;
  mprintf("\tSorting the following sets:\n");
  setsToSort.List();
# ifdef MPI
  // Number of sets to sort should be equal to # members I am responsible for.
  if (Parallel::N_Ens_Members() != (int)setsToSort.size()) {
    rprinterr("Internal Error: Number of ensemble members (%i) != # sets to sort (%zu)\n",
               Parallel::N_Ens_Members(), setsToSort.size());
    return 1;
  }
# endif

  DataSet::DataType dtype = setsToSort[0]->Type();
  unsigned int maxSize = 0;
  for (DataSetList::const_iterator ds = setsToSort.begin(); ds != setsToSort.end(); ++ds) {
    if ((*ds)->Size() < 1) { //TODO check sizes match
      rprinterr("Error: Set '%s' is empty.\n", (*ds)->legend());
      err = 1;
      break;
    }
    if (ds == setsToSort.begin())
      maxSize = (*ds)->Size();
    else if ((*ds)->Size() < maxSize) {
      rprintf("Warning: Set '%s' has fewer frames (%zu) than previous set(s) (%u)\n"
              "Warning: Only using the first %zu frames of all sets.\n",
              (*ds)->legend(), (*ds)->Size(), maxSize, (*ds)->Size());
      maxSize = (unsigned int)(*ds)->Size();
    } else if ((*ds)->Size() > maxSize) {
      rprintf("Warning: Set '%s' has more frames (%zu) than previous set(s) (%u)\n"
              "Warning: Only using the first %u frames of all sets.\n",
              (*ds)->legend(), (*ds)->Size(), maxSize, maxSize);
    }
    if (dtype != (*ds)->Type()) {
      rprinterr("Error: Set '%s' has different type than first set.\n", (*ds)->legend());
      err = 1;
      break;
    }
  }
  if (CheckError(err)) return 1;

# ifdef MPI
  unsigned int threadSize = maxSize;
  Parallel::EnsembleComm().AllReduce( &maxSize, &threadSize, 1, MPI_UNSIGNED, MPI_MIN );
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
  if (dtype != DataSet::PH_EXPL) {
    rprinterr("Error: Only works for pH REMD data for now.\n");
    return 1;
  }

  err = Sort_pH_Data( setsToSort, OutputSets, maxSize );

  return err;
}

// Exec_SortEnsembleData::Help()
void Exec_SortEnsembleData::Help() const
{
  mprintf("\t<dset arg0> [<dset arg1> ...]\n"
          "  Sort unsorted data sets. Currently only works for constant pH REMD data.\n");
}


// Exec_SortEnsembleData::Execute()
Exec::RetType Exec_SortEnsembleData::Execute(CpptrajState& State, ArgList& argIn)
{
  debug_ = State.Debug();
  DataSetList setsToSort;
  std::string dsarg = argIn.GetStringNext();
  while (!dsarg.empty()) {
    setsToSort += State.DSL().GetMultipleSets( dsarg );
    dsarg = argIn.GetStringNext();
  }

  int err = 0;
# ifdef MPI
  // For now, require ensemble mode in parallel.
  if (Parallel::EnsembleComm().IsNull()) {
    rprinterr("Error: Data set ensemble sort requires ensemble mode in parallel.\n");
    return CpptrajState::ERR;
  }
  // Only TrajComm masters have complete data.
  if (Parallel::TrajComm().Master()) {
# endif
    DataSetList OutputSets;
    err = SortData( setsToSort, OutputSets );
    if (err == 0) {
      // Remove unsorted sets. 
      for (DataSetList::const_iterator ds = setsToSort.begin(); ds != setsToSort.end(); ++ds)
        State.DSL().RemoveSet( *ds );
      // Add sorted sets.
      for (DataSetList::const_iterator ds = OutputSets.begin(); ds != OutputSets.end(); ++ds)
        State.DSL().AddSet( *ds );
      // Since sorted sets have been transferred to master DSL, OutputSets now
      // just has copies.
      OutputSets.SetHasCopies( true );
      mprintf("\tSorted sets:\n");
      OutputSets.List();
    }
# ifdef MPI
  }
  if (Parallel::World().CheckError( err ))
# else
  if (err != 0) 
# endif
    return CpptrajState::ERR;
  return CpptrajState::OK;
}
