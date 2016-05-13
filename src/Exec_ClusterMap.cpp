#include "Exec_ClusterMap.h"
#include "CpptrajStdio.h"
#include "ClusterMap.h"
#include "DataSet_MatrixFlt.h"

Exec_ClusterMap::Exec_ClusterMap() : Exec(HIDDEN) {}

void Exec_ClusterMap::Help() const {
  mprintf("\t<2D set> [minpoints <#>] [epsilon <epsilon>] [name <setname>]\n"
          "\t[out <outfile>] [cmapdetail]\n"
          "  Cluster regions of the given 2D data set.\n");
}

// Exec_ClusterMap::Execute()
Exec::RetType Exec_ClusterMap::Execute(CpptrajState& State, ArgList& argIn)
{
  bool cmap_square = !argIn.hasKey("cmapdetail");
  int minPoints = argIn.getKeyInt("minpoints", 4);
  double epsilon = argIn.getKeyDouble("epsilon", 10.0);
  ClusterMap CMAP;
  if (CMAP.Init(epsilon, minPoints)) return CpptrajState::ERR;
  mprintf("\tminpoints= %i, epsilon= %f\n", CMAP.MinPoints(), CMAP.Epsilon());

  std::string dsname = argIn.GetStringKey("name");
  DataFile* outfile = State.DFL().AddDataFile( argIn.GetStringKey("out"), argIn );
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

  if (CMAP.DoCluster( matrix )) return CpptrajState::ERR;

  mprintf("\t%zu clusters:\n", CMAP.Clusters().size());

  // Process clusters.
  Dimension const& ColDim = matrix.Dim(0);
  Dimension const& RowDim = matrix.Dim(1);
  for (ClusterMap::const_iterator CL = CMAP.Clusters().begin(); CL != CMAP.Clusters().end(); ++CL)
  {
    // Write cluster map
    if (cmap_square) {
      int MAXR = CL->MaxRow() + 1; // Want up to and including max
      int MAXC = CL->MaxCol() + 1;
      for (int row = CL->MinRow(); row != MAXR; row++)
        for (int col = CL->MinCol(); col != MAXC; col++)
          output.SetElement( col, row, CL->Cnum() );
    } else {
      for (ClusterMap::Iarray::const_iterator pt = CL->Points().begin();
                                              pt != CL->Points().end(); ++pt)
        output[*pt] = CL->Cnum();
    }
    mprintf("\t %6i: %8zu points, Rows %6g - %6g  Cols %6g - %6g  Avg= %g\n",
            CL->Cnum(), CL->Points().size(),
            RowDim.Coord(CL->MinRow()), RowDim.Coord(CL->MaxRow()),
            ColDim.Coord(CL->MinCol()), ColDim.Coord(CL->MaxCol()), CL->Avg());
  }

  return CpptrajState::OK;
}
