#include <algorithm> //sort, unique
#include <cmath>
#include "Exec_ClusterMap.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
/*
#ifdef _OPENMP
#  include <omp.h>
#endif
*/

Exec_ClusterMap::Exec_ClusterMap() : Exec(HIDDEN),
  epsilon_(0.0),
  epsilon2_(0.0),
  Avg_(0.0),
  minPoints_(0),
  nClusters_(0)
{}

void Exec_ClusterMap::Help() const {
  mprintf("\t<2D set> [minpoints <#>] [epsilon <epsilon>] [name <setname>]\n"
          "\t[out <outfile>] [cmapdetail]\n"
          "  Cluster regions of the given 2D data set.\n");
}

static inline void IdxToColRow(int idx, int ncols, int& col, int& row) {
  col = idx % ncols;
  row = idx / ncols;
}

// Exec_ClusterMap::Execute()
Exec::RetType Exec_ClusterMap::Execute(CpptrajState& State, ArgList& argIn)
{
  cmap_square_ = !argIn.hasKey("cmapdetail");
  minPoints_ = argIn.getKeyInt("minpoints", 4);
  epsilon_ = argIn.getKeyDouble("epsilon", 10.0);
  epsilon2_ = epsilon_ * epsilon_;
  // Based on epsilon, determine max # rows/cols we will have to go. Round up.
  idx_offset_ = (int)ceil(epsilon_);
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
/*
# ifdef _OPENMP
  numthreads_ = 0;
# pragma omp parallel
  {
  if (omp_get_thread_num() == 0)
    numthreads_ = omp_get_num_threads();
  }
  mprintf("\tParallelizing calc with %i threads.\n", numthreads_);
  thread_neighbors_.resize( numthreads_ );
# endif
*/
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
  int maxCol, maxRow;
  IdxToColRow( maxIdx, matrix.Ncols(), maxCol, maxRow );
  mprintf("\t%zu elements, max= %f at index %u (%i, %i), Avg= %f\n",
          matrix.Size(), maxVal, maxIdx, maxCol, maxRow, Avg_);
  mprintf("\tPoints below %f will be treated as noise.\n", Avg_);

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
  int iterations = 0;
# ifdef TIMER
  t_overall_.Start();
# endif
  for (int point = 0; point != (int)matrix.Size(); point++)
  {
    if (!Visited[point])
    {
      progress.Update(iterations++);
      // Mark this point as visited.
      Visited[point] = true;
      double val = matrix.GetElement(point);
      // If this point is less than the cutoff just mark it as noise.
      if (val < Avg_)
        Status[point] = NOISE;
      else
      {
        // Determine how many other points are near this point.
#       ifdef TIMER
        t_query1_.Start();
#       endif
        RegionQuery( NeighborPts, val, point, matrix );
#       ifdef TIMER
        t_query1_.Stop();
#       endif
#       ifdef DEBUG_CLUSTERMAP
        mprintf("\tPoint %i has %u neighbors:", point, NeighborPts.size());
#       endif
        // If # of neighbors less than cutoff, noise; otherwise cluster.
        if ((int)NeighborPts.size() < minPoints_)
        {
#         ifdef DEBUG_CLUSTERMAP
          mprintf(" NOISE\n");
#         endif
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
              progress.Update(iterations++);
#             ifdef DEBUG_CLUSTERMAP
              //mprintf(" %i", neighbor_pt + 1);
#             endif
              // Mark this neighbor as visited
              Visited[neighbor_pt] = true;
              // Determine how many other points are near this neighbor
#             ifdef TIMER
              t_query2_.Start();
#             endif
              RegionQuery( Npts2, neighbor_val, neighbor_pt, matrix );
#             ifdef TIMER
              t_query2_.Stop();
#             endif
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
#         ifdef DEBUG_CLUSTERMAP
          mprintf(" %zu frames,", cluster_frames.size());
#         endif
          std::sort(cluster_frames.begin(), cluster_frames.end());
          Iarray::iterator it = std::unique(cluster_frames.begin(),
                                            cluster_frames.end());
          cluster_frames.resize( std::distance(cluster_frames.begin(),it) );
          // Add cluster to the list
#         ifdef DEBUG_CLUSTERMAP
          mprintf(" %zu unique frames.\n", cluster_frames.size());
#         endif
          AddCluster( cluster_frames, matrix );
        }
      } // END value > cutoff
    } // END if not visited
  } // END loop over matrix points
# ifdef TIMER
  t_overall_.Stop();
  t_query1_.WriteTiming(2, "Region Query 1:", t_overall_.Total());
  t_query2_.WriteTiming(2, "Region Query 2:", t_overall_.Total());
  t_overall_.WriteTiming(1, "Overall:");
# endif
  mprintf("\t%zu clusters:\n", clusters_.size());

  // Sort by number of points
  std::sort(clusters_.begin(), clusters_.end());
  Dimension const& ColDim = matrix.Dim(0);
  Dimension const& RowDim = matrix.Dim(1);
  // Renumber clusters.
  int cnum = 0;
  for (Carray::iterator CL = clusters_.begin(); CL != clusters_.end(); ++CL)
  {
    CL->SetCnum( cnum );
    // Write cluster map
    if (cmap_square_) {
      int MAXR = CL->MaxRow() + 1; // Want up to and including max
      int MAXC = CL->MaxCol() + 1;
      for (int row = CL->MinRow(); row != MAXR; row++)
        for (int col = CL->MinCol(); col != MAXC; col++)
          output.SetElement( col, row, cnum );
    } else {
      for (Iarray::const_iterator pt = CL->Points().begin(); pt != CL->Points().end(); ++pt)
        output[*pt] = cnum;
    }
    mprintf("\t %6i: %8zu points, Rows %6g - %6g  Cols %6g - %6g  Avg= %g\n",
            CL->Cnum(), CL->Points().size(),
            RowDim.Coord(CL->MinRow()), RowDim.Coord(CL->MaxRow()),
            ColDim.Coord(CL->MinCol()), ColDim.Coord(CL->MaxCol()), CL->Avg());
    cnum++;
  }

  return CpptrajState::OK;
}

void Exec_ClusterMap::RegionQuery(Iarray& NeighborPts, double val, int point,
                                  DataSet_2D const& matrix)
{
  NeighborPts.clear();
  int ncols = (int)matrix.Ncols();
  int nrows = (int)matrix.Nrows();
  int point_col, point_row;
  IdxToColRow( point, ncols, point_col, point_row );
  int row_beg = std::max(0,     point_row - idx_offset_);
  int row_end = std::min(nrows, point_row + idx_offset_ + 1);
  int col_beg = std::max(0,     point_col - idx_offset_);
  int col_end = std::min(ncols, point_col + idx_offset_ + 1);

  for (int row = row_beg; row != row_end; row++)
  {
    int idx = row * ncols;
    double dr = (double)(point_row - row);
    for (int col = col_beg; col != col_end; col++)
    {
      int otherpoint = idx + col;
      if (point != otherpoint) {
        double other_val = matrix.GetElement( otherpoint );
        if (other_val > Avg_)
        {
          // Distance calculation.
          double dv = val - other_val;
          double dc = (double)(point_col - col);
          double dist2 = dv*dv + dr*dr + dc*dc;
          if ( dist2 < epsilon2_ )
            NeighborPts.push_back( otherpoint );
        }
      }
    }
  }
}

/*
// Exec_ClusterMap::RegionQuery()
void Exec_ClusterMap::RegionQuery(Iarray& NeighborPts, double val, int point,
                                  DataSet_2D const& matrix)
{
  NeighborPts.clear();
# ifdef _OPENMP
  for (std::vector<Iarray>::iterator it = thread_neighbors_.begin();
                                     it != thread_neighbors_.end(); ++it)
    it->clear();
# endif
  int point_col, point_row;
  IdxToColRow( point, matrix.Ncols(), point_col, point_row );
  int msize = (int)matrix.Size();
  int otherpoint;

# ifdef _OPENMP
  int mythread = 0;
# pragma omp parallel private(otherpoint, mythread)
  {
  mythread = omp_get_thread_num();
# pragma omp for
# endif
  for (otherpoint = 0; otherpoint < msize; otherpoint++)
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
#       ifdef _OPENMP
        thread_neighbors_[mythread].push_back( otherpoint );
#       else
        NeighborPts.push_back( otherpoint );
#       endif
    }
  }
# ifdef _OPENMP
  } // END pragma omp parallel
  for (std::vector<Iarray>::const_iterator it = thread_neighbors_.begin();
                                           it != thread_neighbors_.end(); ++it)
    for (Iarray::const_iterator pt = it->begin(); pt != it->end(); ++pt)
      NeighborPts.push_back( *pt );
# endif
}
*/
// Exec_ClusterMap::AddCluster()
void Exec_ClusterMap::AddCluster(Iarray const& points, DataSet_2D const& matrix)
{
# ifdef DEBUG_CLUSTERMAP_ADDCLUSTER
  mprintf("Cluster %i (%zu):", nClusters_, points.size());
# endif
  int ncols = (int)matrix.Ncols();
  int min_col, min_row;
  IdxToColRow(points.front(), ncols, min_col, min_row);
  int max_col = min_col;
  int max_row = min_row;
  int row, col;
  double cavg = 0.0;
  for (Iarray::const_iterator pt = points.begin(); pt != points.end(); ++pt) {
#   ifdef DEBUG_CLUSTERMAP_ADDCLUSTER
    mprintf(" %i", *pt);
#   endif
    // Determine min/max row/col and average of all points
    IdxToColRow( *pt, ncols, col, row );
    min_col = std::min(col, min_col);
    max_col = std::max(col, max_col);
    min_row = std::min(row, min_row);
    max_row = std::max(row, max_row);
    cavg += matrix.GetElement( *pt );
  }
  cavg /= (double)points.size();
  clusters_.push_back( Cluster(points, cavg, nClusters_, min_col, max_col, min_row, max_row) );
# ifdef DEBUG_CLUSTERMAP_ADDCLUSTER
  mprintf("\n");
# endif
  nClusters_++;
}
