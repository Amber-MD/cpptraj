#include <algorithm> // sort
#include <cmath> // ceil
#include "ClusterMap.h"
#include "CpptrajStdio.h"
#include "Constants.h" // SMALL
#include "ProgressBar.h"
#include "ProgressTimer.h"

ClusterMap::ClusterMap() :
  epsilon_(0.0),
  epsilon2_(0.0),
  Avg_(0.0),
  minPoints_(-1),
  idx_offset_(0)
{}

int ClusterMap::Init(double epsilonIn, int minPointsIn)
{
  epsilon_ = epsilonIn;
  if (epsilon_ < Constants::SMALL)
  {
    mprinterr("Error: Epsilon must be > 0.0\n");
    return 1;
  }
  epsilon2_ = epsilon_ * epsilon_;
  // Based on epsilon, determine max # rows/cols we will have to go. Round up.
  idx_offset_ = (int)ceil(epsilon_);

  minPoints_ = minPointsIn;
  if (minPoints_ == 0) {
    mprinterr("Error: Minimum number of points must be > 0\n");
    return 1;
  }
  return 0;
}

static inline void IdxToColRow(int idx, int ncols, int& col, int& row) {
  col = idx % ncols;
  row = idx / ncols;
}

int ClusterMap::DoCluster(DataSet_2D const& matrix)
{
# ifdef TIMER
  t_overall_.Start();
# endif
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

  if ( DoDBSCAN( matrix ) ) return 1;

  mprintf("\t%zu clusters:\n", clusters_.size());

  // Sort by number of points
  std::sort(clusters_.begin(), clusters_.end());

  // Renumber clusters.
  int cnum = 0;
  for (Carray::iterator CL = clusters_.begin(); CL != clusters_.end(); ++CL)
    CL->SetCnum( cnum++ );
# ifdef TIMER
  t_overall_.Stop();
  t_query1_.WriteTiming(2, "Region Query 1:", t_overall_.Total());
  t_query2_.WriteTiming(2, "Region Query 2:", t_overall_.Total());
  t_overall_.WriteTiming(1, "Overall:");
# endif

  return 0;
}

#define UNCLASSIFIED -2
#define NOISE -1

int ClusterMap::DoDBSCAN(DataSet_2D const& matrix)
{
  Status_.assign( matrix.Size(), UNCLASSIFIED );

  // SetOfPoints is UNCLASSIFIED
  int ClusterId = 0;
  ProgressBar progress(matrix.Size());
  ProgressTimer ptimer(matrix.Size());
  for (unsigned int idx = 0; idx != matrix.Size(); idx++)
  {
    progress.Update(idx);
    ptimer.Remaining(idx);
    //Point := SetOfPoints.get(i);
    //IF Point.ClId = UNCLASSIFIED THEN
    if ( Status_[idx] == UNCLASSIFIED )
    {
      if (matrix.GetElement(idx) < Avg_) // NOTE: Not part of original DBSCAN algorithm
        Status_[idx] = NOISE;
      else
      {
        //IF ExpandCluster(SetOfPoints, Point, ClusterId, Eps, MinPts) THEN
        if (ExpandCluster(idx, ClusterId, matrix))
        //ClusterId := nextId(ClusterId)
          ClusterId++;
      }
    }
  }

  mprintf("DEBUG: %i clusters.\n", ClusterId);
  if (ClusterId > 0) {
    std::vector<Iarray> C0( ClusterId );
    for (unsigned int idx = 0; idx != Status_.size(); idx++)
    {
      int point = Status_[idx];
      if (point == UNCLASSIFIED)
        mprintf("Warning: point %u was unclassified.\n", idx);
      else if (point != NOISE)
        C0[ point ].push_back( idx );
    }   
    for (std::vector<Iarray>::const_iterator it = C0.begin(); it != C0.end(); ++it)
      AddCluster( *it, matrix );
  }   
  return 0;
}

bool ClusterMap::ExpandCluster(unsigned int point, int ClusterId, DataSet_2D const& matrix)
{
  //seeds:=SetOfPoints.regionQuery(Point,Eps);
  //Iarray seeds, result;
# ifdef TIMER
  t_query1_.Start();
# endif
  RegionQuery(seeds_, point, matrix);
# ifdef TIMER
  t_query1_.Stop();
# endif

  //IF seeds.size<MinPts THEN // no core point
  if ((int)seeds_.size() < minPoints_)
  {
    //SetOfPoint.changeClId(Point,NOISE);
    Status_[point] = NOISE;
    //RETURN False;
    return false;
  }
  else
  {
    // all points in seeds are density-reachable from Point
    //SetOfPoints.changeClIds(seeds,ClId);
    Status_[point] = ClusterId;
    for (Iarray::const_iterator pt = seeds_.begin(); pt != seeds_.end(); ++pt)
      Status_[*pt] = ClusterId;
    //seeds.delete(Point);
    //WHILE seeds <> Empty DO
    unsigned int endIdx = seeds_.size();
    for (unsigned int idx = 0; idx < endIdx; idx++)
    {
      //currentP := seeds.first();
      int otherpoint = seeds_[idx];
      //result := SetOfPoints.regionQuery(currentP, Eps);
#     ifdef TIMER
      t_query2_.Start();
#     endif
      RegionQuery(result_, otherpoint, matrix);
#     ifdef TIMER
      t_query2_.Stop();
#     endif
      //IF result.size >= MinPts THEN
      if ( (int)result_.size() >= minPoints_ )
      {
        //FOR i FROM 1 TO result.size DO
        //  resultP := result.get(i);
        //  IF resultP.ClId IN {UNCLASSIFIED, NOISE} THEN
        //    IF resultP.ClId = UNCLASSIFIED THEN
        //      seeds.append(resultP);
        //    END IF;
        //    SetOfPoints.changeClId(resultP,ClId);
        //  END IF; // UNCLASSIFIED or NOISE
        //END FOR;
        for (Iarray::const_iterator rt = result_.begin(); rt != result_.end(); ++rt)
        {
          if (Status_[*rt] == UNCLASSIFIED || Status_[*rt] == NOISE)
          {
            if (Status_[*rt] == UNCLASSIFIED)
            {
              seeds_.push_back( *rt );
              endIdx = seeds_.size();
            }
            Status_[*rt] = ClusterId;
          }
        }
      }
      //END IF; // result.size >= MinPts
      //seeds.delete(currentP);
    }
    //END WHILE; // seeds <> Empty
    return true;
  }
}

void ClusterMap::RegionQuery(Iarray& NeighborPts, int point, DataSet_2D const& matrix)
{
  double val = matrix.GetElement(point);
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

void ClusterMap::AddCluster(Iarray const& points, DataSet_2D const& matrix)
{
  int cnum = (int)clusters_.size();
# ifdef DEBUG_CLUSTERMAP_ADDCLUSTER
  mprintf("Cluster %i (%zu):", cnum, points.size());
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
  clusters_.push_back( Cluster(points, cavg, cnum, min_col, max_col, min_row, max_row) );
# ifdef DEBUG_CLUSTERMAP_ADDCLUSTER
  mprintf("\n");
# endif
}

#undef UNCLASSIFIED
#undef NOISE
