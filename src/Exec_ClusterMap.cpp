#include <algorithm> //sort, unique
#include "Exec_ClusterMap.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"

Exec_ClusterMap::Exec_ClusterMap() : Exec(HIDDEN),
  epsilon_(0.0),
  epsilon2_(0.0),
  Avg_(0.0),
  minPoints_(0),
  nClusters_(0)
{}

void Exec_ClusterMap::Help() const {
  mprintf("\t<2d set>\n");
}

static inline void IdxToColRow(int idx, int ncols, int& col, int& row) {
  col = idx % ncols;
  row = idx / ncols;
}

// Exec_ClusterMap::Execute()
Exec::RetType Exec_ClusterMap::Execute(CpptrajState& State, ArgList& argIn)
{
  minPoints_ = argIn.getKeyInt("minPoints", 4);
  epsilon_ = argIn.getKeyDouble("epsilon", 10.0);
  epsilon2_ = epsilon_ * epsilon_;
  std::string dsname = argIn.GetStringKey("name");
  DataFile* outfile = State.DFL().AddDataFile( argIn.GetStringKey("out"), argIn );
  mprintf("\tminpoints= %i, epsilon= %f\n", minPoints_, epsilon_);
  DataSet* dsIn = State.DSL().GetDataSet( argIn.GetStringNext() );
  if (dsIn == 0) return CpptrajState::ERR;
  mprintf("\tSet '%s'\n", dsIn->legend());
  if (dsIn->Group() != DataSet::MATRIX_2D) {
    mprinterr("Error: Set is not 2D.\n");
    return CpptrajState::ERR;
  }
  if (dsIn->Size() < 1) {
    mprinterr("Error: Set is empty.\n");
    return CpptrajState::ERR;
  }
  DataSet_2D const& matrix = static_cast<DataSet_2D const&>( *dsIn );

  // Output set
  if (dsname.empty())
    dsname = State.DSL().GenerateDefaultName("cmap");
  DataSet* dsOut = State.DSL().AddSet( DataSet::MATRIX_FLT, dsname );
  if (dsOut == 0) return CpptrajState::ERR;
  if (outfile != 0) outfile->AddDataSet( dsOut );
  DataSet_MatrixFlt& output = static_cast<DataSet_MatrixFlt&>( *dsOut );
  output.Allocate2D( matrix.Ncols(), matrix.Nrows() );
  std::fill( output.begin(), output.end(), -1.0 );

  unsigned int ncols = matrix.Ncols();
  // Go through the set, calculate the average. Also determine the max.
  double maxVal = matrix.GetElement(0);
  unsigned int maxIdx = 0;
  Avg_ = 0.0;
  for (unsigned int i = 1; i != matrix.Size(); i++) {
    double val = matrix.GetElement(i);
    if (val > maxVal) {
      maxVal = val;
      maxIdx = i;
    }
    Avg_ += val;
  }
  Avg_ /= (double)matrix.Size();
  unsigned int maxCol = maxIdx % ncols;
  unsigned int maxRow = maxIdx / ncols;
  mprintf("\t%zu elements, max= %f at index %u (%u, %u), Avg= %f\n",
          matrix.Size(), maxVal, maxIdx, maxCol, maxRow, Avg_);

  // Use DBSCAN-style algorithm to cluster points. Any point less than the
  // average will be considered noise.
  std::vector<bool> Visited( matrix.Size(), false );
  const char UNASSIGNED = 'U';
  const char NOISE = 'N';
  const char INCLUSTER = 'C';
  std::vector<char> Status( matrix.Size(), UNASSIGNED );
  // Loop over all points in matrix
  Iarray NeighborPts;    // Will hold all neighbors of the current point
  Iarray Npts2;          // Will hold neighbors of a neighbor
  Iarray cluster_frames; // Hold indices of current cluster
  ProgressBar progress(matrix.Size());
  for (int point = 0; point != (int)matrix.Size(); point++)
  {
    progress.Update(point);
    if (!Visited[point])
    {
      // Mark this point as visited.
      Visited[point] = true;
      double val = matrix.GetElement(point);
      // If this point is less than the cutoff just mark it as noise.
      if (val < Avg_)
        Status[point] = NOISE;
      else
      {
        // Determine how many other points are near this point.
        RegionQuery( NeighborPts, val, point, matrix );
        mprintf("\tPoint %u\n", point);
        mprintf("\t\t%u neighbors:\n", NeighborPts.size());
        // If # of neighbors less than cutoff, noise; otherwise cluster.
        if ((int)NeighborPts.size() < minPoints_)
        {
          mprintf(" NOISE\n");
          Status[point] = NOISE;
        }
        else
        {
          // Expand cluster
          cluster_frames.clear();
          cluster_frames.push_back( point );
          // NOTE: Use index instead of iterator since NeighborPts may be
          //       modified inside this loop.
          unsigned int endidx = NeighborPts.size();
          for (unsigned int idx = 0; idx < endidx; ++idx)
          {
            int neighbor_pt = NeighborPts[idx];
            double neighbor_val = matrix.GetElement(neighbor_pt);
            if (!Visited[neighbor_pt])
            {
              mprintf(" %i", neighbor_pt + 1);
              // Mark this neighbor as visited
              Visited[neighbor_pt] = true;
              // Determine how many other points are near this neighbor
              RegionQuery( Npts2, neighbor_val, neighbor_pt, matrix );
              if ((int)Npts2.size() >= minPoints_)
              {
                // Add other points to current neighbor list
                NeighborPts.insert( NeighborPts.end(), Npts2.begin(), Npts2.end() );
                endidx = NeighborPts.size();
              }
            }
            // If neighbor is not yet part of a cluster, add it to this one.
            if (Status[neighbor_pt] != INCLUSTER)
            {
              cluster_frames.push_back( neighbor_pt );
              Status[neighbor_pt] = INCLUSTER;
            }
          }
          // Remove duplicate frames
          // TODO: Take care of this in Renumber?
          std::sort(cluster_frames.begin(), cluster_frames.end());
          Iarray::iterator it = std::unique(cluster_frames.begin(),
                                            cluster_frames.end());
          cluster_frames.resize( std::distance(cluster_frames.begin(),it) );
          // Add cluster to the list
          mprintf("\n");
          AddCluster( cluster_frames, output );
        }
      } // END value > cutoff
    } // END if not visited
  } // END loop over matrix points

  return CpptrajState::OK;
}

// Exec_ClusterMap::RegionQuery()
void Exec_ClusterMap::RegionQuery(Iarray& NeighborPts, double val, int point,
                                  DataSet_2D const& matrix)
{
  NeighborPts.clear();
  int point_col, point_row;
  IdxToColRow( point, matrix.Ncols(), point_col, point_row );
  
  for (int otherpoint = 0; otherpoint != (int)matrix.Size(); otherpoint++)
  {
    if (point != otherpoint) {
      // Distance calculation.
      int other_col, other_row;
      IdxToColRow( otherpoint, matrix.Ncols(), other_col, other_row );
      double other_val = matrix.GetElement( otherpoint );
      double dv = val - other_val;
      double dr = (double)(point_row - other_row);
      double dc = (double)(point_col - other_col);
      double dist2 = dv*dv + dr*dr + dc*dc;
      
      if ( dist2 < epsilon2_ )
        NeighborPts.push_back( otherpoint );
    }
  }
}

void Exec_ClusterMap::AddCluster(Iarray const& points, DataSet_MatrixFlt& output) {
  mprintf("Cluster %i:", nClusters_);
  for (Iarray::const_iterator it = points.begin(); it != points.end(); ++it) {
    mprintf(" %i", *it);
    output[*it] = nClusters_;
  }
  mprintf("\n");
  nClusters_++;
}
