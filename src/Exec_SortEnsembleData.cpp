#include "Exec_SortEnsembleData.h"
#include "CpptrajStdio.h"
#include "DataSet_PHREMD_Explicit.h"
#include "DataSet_PHREMD_Implicit.h"
#include "DataSet_pH.h"
#include "StringRoutines.h" // doubleToString

//  Exec_SortEnsembleData::Sort_pH_Data()
int Exec_SortEnsembleData::Sort_pH_Data(DataSetList const& setsToSort, DataSetList& OutputSets,
                                        unsigned int maxFrames)
const
{
  // Cast sets back to DataSet_PHREMD
  typedef std::vector<DataSet_PHREMD*> Parray;
  Parray PHsets;
  for (DataSetList::const_iterator ds = setsToSort.begin(); ds != setsToSort.end(); ++ds)
    PHsets.push_back( (DataSet_PHREMD*)*ds );

  // Gather initial pH data values, ensure no duplicates
  typedef std::vector<double> Darray;
  Darray pHvalues;
# ifdef MPI
  pHvalues.resize( Parallel::Ensemble_Size() );
  Darray phtmp;
  for (Parray::const_iterator ds = PHsets.begin(); ds != PHsets.end(); ++ds)
    phtmp.push_back( (*ds)->Initial_pH() );
  if (comm_.AllGather(&phtmp[0], phtmp.size(), MPI_DOUBLE, &pHvalues[0])) {
    rprinterr("Error: Gathering pH values.\n");
    return 1;
  }
# else
  for (Parray::const_iterator ds = PHsets.begin(); ds != PHsets.end(); ++ds)
    pHvalues.push_back( (*ds)->Initial_pH() );
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
  DataSet_PHREMD::Rarray const& Residues = PHsets[0]->Residues();
  int defaultState = 0;
# ifdef MPI
  if ( PHsets[0]->Type() == DataSet::PH_IMPL)
    defaultState = -1;
# endif
  if (debug_ > 0)
    rprintf("DEBUG: Sorting %u frames for %zu sets, %zu pH values.\n",
            maxFrames, PHsets.size(), sortedPH.size());
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
      out->Resize( maxFrames, defaultState );
    }
  }

  // ---------------------------------------------
  if ( PHsets[0]->Type() == DataSet::PH_EXPL) {
    // Loop over unsorted sets
    for (Parray::const_iterator ds = PHsets.begin(); ds != PHsets.end(); ++ds)
    {
      DataSet_PHREMD_Explicit* in = (DataSet_PHREMD_Explicit*)*ds;
      unsigned int phidx = 0;
      for (unsigned int n = 0; n < maxFrames; n++)
      {
        float phval = in->pH_Values()[n];
        int setidx = pH_map.FindIndex( phval ) * Residues.size();
        //rprintf("DEBUG: %6u Set %10s pH= %6.2f going to %2i\n", n+1, in->legend(), phval, idx);
        //mflush();
        for (unsigned int res = 0; res < in->Residues().size(); res++, setidx++, phidx++)
        {
          DataSet_pH* out = (DataSet_pH*)OutputSets[setidx];
          //if (res == 0 && idx == 0) {
          //  rprintf("DEBUG: Frame %3u res %2u State %2i pH %6.2f\n", 
          //          n, res, in->Res(res).State(n), phval);
          //  mflush();
          //}
          out->SetState(n, in->ResStates()[phidx], in->RecordType(n));
        }
      }
    } // END loop over unsorted sets
#   ifdef MPI
    // Now we need to reduce down each set onto the thread where it belongs.
    if (Parallel::World().Size() > 1) {
      for (int idx = 0; idx != (int)OutputSets.size(); idx++) {
        DataSet_pH* out = (DataSet_pH*)OutputSets[idx];
        int ensembleRank = Parallel::MemberEnsCommRank( out->Meta().EnsembleNum() );
        //rprintf("DEBUG: Reduce set %s to rank %i\n", out->legend(), ensembleRank);
        out->Reduce( comm_, ensembleRank );
      }
      // Remove sets that do not belong on this rank
      for (int idx = (int)OutputSets.size() - 1; idx > -1; idx--) {
        DataSet* out = OutputSets[idx];
        int ensembleRank = Parallel::MemberEnsCommRank( out->Meta().EnsembleNum() );
        if (ensembleRank != comm_.Rank()) {
          //rprintf("DEBUG: Remove set %s (%i) from rank %i\n", out->legend(),
          //        idx, comm_.Rank());
          OutputSets.RemoveSet( out );
        }
      }
    }
#   endif
  // ---------------------------------------------
  } else if ( PHsets[0]->Type() == DataSet::PH_IMPL) {
#   ifdef MPI
    typedef std::vector<int> Iarray;
    typedef std::vector<Iarray> Iarray2;
    typedef std::vector<bool> Barray;
    // True if I have this pH value
    Barray isMyPh( sortedPH.size(), false );
    // Which rank in ensemble has which pH
    Iarray pHrank( sortedPH.size(), 0 );
    for (int phidx = 0; phidx != (int)sortedPH.size(); phidx++)
    {
      int ensembleRank = Parallel::MemberEnsCommRank( phidx );
      pHrank[phidx] = ensembleRank;
      isMyPh[phidx] = (ensembleRank == comm_.Rank());
    }
    // DEBUG
    for (unsigned int idx = 0; idx != pHrank.size(); idx++)
      mprintf("\tpH %6.2f on rank %i\n", sortedPH[idx], pHrank[idx]);
    // Hold frame-residue-state for each pH not on this rank.
    Iarray2 FrmResState( sortedPH.size() );
    // Loop over unsorted sets
    for (Parray::const_iterator ds = PHsets.begin(); ds != PHsets.end(); ++ds)
    {
      DataSet_PHREMD_Implicit* in = (DataSet_PHREMD_Implicit*)*ds;
      // Loop over frames
      for (unsigned int n = 0; n < maxFrames; n++)
      {
        DataSet_PHREMD_Implicit::Record const& Rec = in->Records()[n];
        float phval = Rec.pH();
        int phidx = pH_map.FindIndex( phval );
        if (isMyPh[phidx]) {
          // This pH belongs to me. Set value.
          int setidx = phidx * Residues.size();
          if (Rec.RecType() == Cph::FULL_RECORD) {
            // Info for all residues
            for (unsigned int res = 0; res < in->Residues().size(); res++, setidx++)
              ((DataSet_pH*)OutputSets[setidx])->SetState(n, Rec.ResStates()[res], Rec.RecType());
          } else if (Rec.RecType() > -1) {
            // Info for single residue, record type for all residues
            //rprintf("\tSetting my pH %6.2f frame %8i state %2i idx %6i res %6i '%s'\n", sortedPH[phidx], n, Rec.ResStates()[0], setidx, Rec.RecType(), OutputSets[setidx+Rec.RecType()]->legend());
            for (int res = 0; res < (int)in->Residues().size(); res++, setidx++)
              if (res == Rec.RecType())
                ((DataSet_pH*)OutputSets[setidx])->SetState(n, Rec.ResStates()[0], Rec.RecType());
              else
                ((DataSet_pH*)OutputSets[setidx])->SetRecType(n, Rec.RecType());
          }
        } else {
          // This pH belongs to another rank. Save it.
          if (Rec.RecType() > -1) {
            // Info for a single residue present
            FrmResState[phidx].push_back( n );
            FrmResState[phidx].push_back( Rec.RecType() );
            FrmResState[phidx].push_back( Rec.ResStates()[0] );
          } else {
            // Info for all residues present
            FrmResState[phidx].push_back( n );
            FrmResState[phidx].push_back( Rec.RecType() );
            for (unsigned int res = 0; res < in->Residues().size(); res++)
              FrmResState[phidx].push_back( Rec.ResStates()[res] );
          }
        }
      } // END loop over frames
    } // END loop over sets
    // DEBUG
/*
    comm_.Barrier();
    for (int rank = 0; rank < comm_.Size(); rank++)
    {
      if (rank == comm_.Rank())
      {
        for (unsigned int phidx = 0; phidx != sortedPH.size(); phidx++) {
          rprintf("DEBUG: pH %6.2f: %8s %6s %2s\n", sortedPH[phidx], "Frm", "Res", "St");
          Iarray const& FRS = FrmResState[phidx];
          unsigned int idx = 0;
          while (idx < FRS.size()) {
            int rec = FRS[idx+1];
            if (rec > -1) {
              rprintf("                  %8i %6i %2i\n", FRS[idx], rec, FRS[idx+2]);
              idx += 3;
            } else {
              rprintf("                  %8i %6i All Residues\n", FRS[idx], rec);
              idx += (2 + Residues.size());
            }
          }
        }
      }
      comm_.Barrier();
    }
*/
    // Communicate states to other ranks
    typedef std::vector<unsigned int> Uarray;
    Uarray sizeOnRank( comm_.Size() );
    for (unsigned int phidx = 0; phidx != sortedPH.size(); phidx++)
    {
      // Each rank says how many frames of this pH they have collected and
      // send to rank the pH belongs to.
      unsigned int nph = FrmResState[phidx].size();
      comm_.Gather(&nph, 1, MPI_UNSIGNED, &sizeOnRank[0], pHrank[phidx]);
      if (pHrank[phidx] == comm_.Rank())
      {
        // This pH belongs to me. I should have no frames at this pH.
        if (sizeOnRank[comm_.Rank()] > 0) {
          rprinterr("Internal Error: Rank has frames to communicate at its pH\n");
          Parallel::Abort(1);
        }
        unsigned int totalSize = 0;
        for (unsigned int idx = 0; idx != sizeOnRank.size(); idx++) {
          totalSize += sizeOnRank[idx];
          //rprintf("DEBUG: Rank %4u has %8u frames of pH %6.2f\n",
          //        idx, sizeOnRank[idx], sortedPH[phidx]);
        }
        //rprintf("DEBUG: Total incoming size: %u\n", totalSize);
        FrmResState[phidx].resize( totalSize );
        // Receive state info for this pH from other ranks
        int* frsArray = &(FrmResState[phidx][0]);
        for (int rank = 0; rank != comm_.Size(); rank++) {
          if (rank != comm_.Rank()) {
            comm_.Recv(frsArray, sizeOnRank[rank], MPI_INT, rank, 1600+rank);
            frsArray += sizeOnRank[rank];
          }
        }
      } else {
        // This pH belongs to another rank. Send my info there.
        int* frsArray = &(FrmResState[phidx][0]);
        comm_.Send(frsArray, nph, MPI_INT, pHrank[phidx], 1600+comm_.Rank());
      }
      comm_.Barrier();
    }
    // Fill in state info
    std::vector<DataSet*> ToRemove;
    for (unsigned int phidx = 0; phidx != sortedPH.size(); phidx++)
    {
      int setidx = phidx * Residues.size();
      if (pHrank[phidx] == comm_.Rank())
      {
        Iarray const& FRS = FrmResState[phidx];
        // This pH belongs to me. Fill in the information received from
        // other ranks.
        unsigned int idx = 0;
        while (idx < FRS.size()) {
          int rec = FRS[idx+1];
          if (rec > -1) {
            // Info for single residue, record type for all residues
            //rprintf("\tSetting pH %6.2f frame %8i state %2i idx %6i res %6i '%s'\n", sortedPH[phidx], FRS[idx], FRS[idx+2], setidx, rec, OutputSets[setidx+rec]->legend());
            for (int res = 0; res != (int)Residues.size(); res++) {
              DataSet_pH* out = (DataSet_pH*)OutputSets[setidx + res];
              if (rec == res)
                out->SetState( FRS[idx], FRS[idx+2], rec );
              else
                out->SetRecType( FRS[idx], rec );
            }
            idx += 3;
          } else {
            //rprintf("                  %8i %6i All Residues\n", FRS[idx], rec);
            int frm = FRS[idx];
            idx += 2;
            for (unsigned int res = 0; res != Residues.size(); res++, idx++) {
              DataSet_pH* out = (DataSet_pH*)OutputSets[setidx + res];
              out->SetState( frm, FRS[idx], rec );
            }
          }
        }
        // Fill in any remaining data. FIXME safe to assume first frame is set?
        for (unsigned int res = 0; res != Residues.size(); res++) {
          DataSet_pH* out = (DataSet_pH*)OutputSets[setidx + res];
          for (unsigned int n = 1; n < maxFrames; n++)
            if (out->State(n) == -1)
              out->SetState(n, out->State(n-1), out->RecordType(n));
        }
      } else {
        // This pH does not belong to me. Mark associated data sets to be removed.
        for (unsigned int res = 0; res != Residues.size(); res++)
          ToRemove.push_back( OutputSets[setidx + res] );
      }
    }
    // Remove data sets that do not belong to me.
    for (std::vector<DataSet*>::reverse_iterator it = ToRemove.rbegin();
                                                 it != ToRemove.rend(); ++it)
    {
      //rprintf("DEBUG: '%s' does not belong to me.\n", (*it)->legend());
      OutputSets.RemoveSet( *it );
    }
#   else /* if not MPI */
    // Loop over frames
    for (unsigned int n = 0; n < maxFrames; n++)
    {
      // Loop over unsorted sets
      for (Parray::const_iterator ds = PHsets.begin(); ds != PHsets.end(); ++ds)
      {
        DataSet_PHREMD_Implicit* in = (DataSet_PHREMD_Implicit*)*ds;
        DataSet_PHREMD_Implicit::Record const& Rec = in->Records()[n];
        float phval = Rec.pH();
        int setidx = pH_map.FindIndex( phval ) * Residues.size();
        if (Rec.RecType() == Cph::FULL_RECORD) {
          for (unsigned int res = 0; res < in->Residues().size(); res++, setidx++)
          {
            DataSet_pH* out = (DataSet_pH*)OutputSets[setidx];
            //if (res == 0 && idx == 0) {
            //  rprintf("DEBUG: Frame %3u res %2u State %2i pH %6.2f\n", 
            //          n, res, in->Res(res).State(n), phval);
            //  mflush();
            //}
            out->SetState(n, Rec.ResStates()[res], Rec.RecType());
          }
        } else {
          for (int res = 0; res < (int)in->Residues().size(); res++, setidx++)
          {
            DataSet_pH* out = (DataSet_pH*)OutputSets[setidx];
            if (res == Rec.RecType())
              out->SetState(n, Rec.ResStates()[0], Rec.RecType());
            else
              // State for this residue not recorded - use previous state.
              // Should be fine since first state always has all residues.
              out->SetState(n, out->State(n-1), Rec.RecType());
          }
        }
      } // END loop over unsorted sets
    } // END loop over frames
#   endif /* MPI */
  // ---------------------------------------------
  } else {
    return 1; // Sanity check
  }
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
  comm_.AllReduce( &maxSize, &threadSize, 1, MPI_UNSIGNED, MPI_MIN );
  typedef std::vector<int> Iarray;
  Iarray Dtypes( comm_.Size(), -1 );
  if ( comm_.AllGather( &dtype, 1, MPI_INT, &Dtypes[0] ) ) return 1;
  for (int rank = 1; rank < comm_.Size(); rank++)
    if (Dtypes[0] != Dtypes[rank]) {
      rprinterr("Error: Set types on rank %i do not match types on rank 0.\n", rank);
      err = 1;
      break;
    }
  if (comm_.CheckError( err )) return 1;
# endif

  // Only work for pH data for now.
  if (dtype != DataSet::PH_EXPL && dtype != DataSet::PH_IMPL) {
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
  if (!Parallel::EnsembleIsSetup()) {
    rprinterr("Error: Data set ensemble sort requires ensemble mode in parallel.\n");
    return CpptrajState::ERR;
  }
  // Only TrajComm masters have complete data.
  if (Parallel::TrajComm().Master()) {
    comm_ = Parallel::MasterComm();
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
