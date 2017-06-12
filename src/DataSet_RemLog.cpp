#include "DataSet_RemLog.h"
#include "CpptrajStdio.h"

DataSet_RemLog::DataSet_RemLog() :
  // 0 dim indicates DataSet-specific write 
  DataSet(REMLOG, GENERIC, TextFormat(TextFormat::DOUBLE, 10, 4), 0)
{}

/** Setup a single 1-dimension group.
  * \param group_size Expected number of replicas.
  */
void DataSet_RemLog::SetupDim1Group(int group_size) {
  groupDims_.clear();
  groupDims_.resize( 1 );
  groupDims_[0].resize( 1 );
  for (int replica = 0; replica < group_size; replica++) {
    int me = replica + 1;
    int l_partner = me - 1;
    if (l_partner < 1) l_partner = group_size;
    int r_partner = me + 1;
    if (r_partner > group_size) r_partner = 1;
    groupDims_[0][0].push_back( GroupReplica(l_partner, me, r_partner) );
  }
}

/** Allocate for 1D REMD. */ // TODO Pass in dimension type, not array
void DataSet_RemLog::AllocateReplicas(int n_replicas, ReplicaDimArray const& repDimIn,
                                       int debugIn)
{
  AllocateReplicas(n_replicas, GdimArray(), repDimIn, debugIn);
}

// DataSet_RemLog::AllocateReplicas()
/** \param n_replicas Total number of replicas across all dimensions.
  * \param gdimIn Array describing layout of replicas; if empty assume simple 1D layout.
  * \param repDimIn Array describing each replica dimension.
  * \param debugIn Debug level; higher means more info printed.
  */
void DataSet_RemLog::AllocateReplicas(int n_replicas, GdimArray const& gdimIn,
                                      ReplicaDimArray const& repDimIn, int debugIn)
{
  ensemble_.clear();
  ensemble_.resize( n_replicas );

  if (gdimIn.empty())
    // Assume 1D
    SetupDim1Group( n_replicas );
  else
    groupDims_ = gdimIn;
  // DEBUG: Print out dimension layout
  if (debugIn > 0) {
    for (GdimArray::const_iterator Dim = groupDims_.begin();
                                   Dim != groupDims_.end(); ++Dim)
    {
      mprintf("Dimension %u:\n", Dim - groupDims_.begin());
      for (GroupDimType::const_iterator Group = Dim->begin();
                                        Group != Dim->end(); ++Group)
      {
        mprintf("\tGroup %u:\n", Group - Dim->begin());
        for (GroupArray::const_iterator Rep = Group->begin();
                                        Rep != Group->end(); ++Rep)
          mprintf("\t\tReplica[%u]= %i (l=%i, r=%i)\n", Rep - Group->begin(), Rep->Me(),
                  Rep->L_partner(), Rep->R_partner());
      }
    }
  }

  // Set up dimension topological info for each replica.
  repInfo_.clear();
  repInfo_.resize( n_replicas ); // [rep][dim]
  for (unsigned int dim = 0; dim != groupDims_.size(); dim++) {
    for (unsigned int grp = 0; grp != groupDims_[dim].size(); grp++) {
      unsigned int topidx = groupDims_[dim][grp].size() - 1;
      for (unsigned int idx = 0; idx != groupDims_[dim][grp].size(); idx++)
      {
        LocationType loc;
        if (idx == 0)
          loc = BOTTOM;
        else if (idx == topidx)
          loc = TOP;
        else
          loc = MIDDLE;
        GroupReplica const& Rep = static_cast<GroupReplica const&>( groupDims_[dim][grp][idx] );
        repInfo_[ Rep.Me() - 1 ].push_back(RepInfo(grp, Rep.L_partner()-1, Rep.R_partner()-1, loc));
      }
    }
  }
  if (debugIn > 0) {
    const char* LOCSTRING[] = {"BOT", "MID", "TOP"};
    for (unsigned int rep = 0; rep != repInfo_.size(); rep++) {
      mprintf("\tReplica %u:", rep);
      for (unsigned int dim = 0; dim != repInfo_[rep].size(); dim++) {
        RepInfo const& RI = static_cast<RepInfo const&>(repInfo_[rep][dim]);
        mprintf(" Dim%u[G=%i L=%i R=%i Loc=%s]", dim, RI.GroupID(), RI.LeftID(), RI.RightID(),
                LOCSTRING[RI.Location()]);
      }
      mprintf("\n");
    }
  }

  repDims_ = repDimIn;
}

/** \return Total number of exchanges based on first replica. */
int DataSet_RemLog::NumExchange() const {
  if (ensemble_.empty())
    return 0;
  else // Each member of the ensemble should have same # exchanges.
    return (int)ensemble_[0].size();
}

/** \return true if all replicas have same number of exchanges. */
bool DataSet_RemLog::ValidEnsemble() const {
  ReplicaEnsemble::const_iterator member = ensemble_.begin();
  size_t first_size = (*member).size();
  for (; member != ensemble_.end(); ++member) {
    if ((*member).size() != first_size) {
      mprinterr("Error: In remlog data set %s size of ensemble member %zu (%zu) !="
                " size of first member (%zu)\n", Meta().Name().c_str(), // TODO: Change to legend
                member - ensemble_.begin() + 1, (*member).size(), first_size);
      return false;
    }
  }
  return true;
}

/** Ensure all replicas have same number of exchanges by setting all to
  * the minimum number of exchanges among all current replicas.
  */
void DataSet_RemLog::TrimLastExchange() {
  if (ensemble_.empty()) return;
  ReplicaEnsemble::iterator member = ensemble_.begin();
  size_t min_size = member->size();
  ++member;
  for (; member != ensemble_.end(); ++member) {
    if (member->size() < min_size) min_size = member->size();
  }
  // Resize all member arrays to minimum
  for (member = ensemble_.begin(); member != ensemble_.end(); ++member)
    member->resize( min_size );
}

// DataSet_RemLog::PrintReplicaStats()
void DataSet_RemLog::PrintReplicaStats() const {
  mprintf("Replica Stats:\n"
          "%-10s %2s %6s %6s %6s %12s %12s %12s S\n", "#Exchange", "#D",
          "RepIdx", "PrtIdx", "CrdIdx", "Temp0", "PE_X1", "PE_X2");
  for (int exchg = 0; exchg < NumExchange(); exchg++) {
    for (int replica = 0; replica < (int)Size(); replica++) {
      ReplicaFrame const& frm = RepFrame(exchg, replica); 
      mprintf("%10u %2i %6i %6i %6i %12.4f %12.4f %12.4f %1i\n", exchg + 1, frm.Dim(),
              frm.ReplicaIdx(), frm.PartnerIdx(), frm.CoordsIdx(), frm.Temp0(), 
              frm.PE_X1(), frm.PE_X2(), (int)frm.Success());
    }
  }
}
