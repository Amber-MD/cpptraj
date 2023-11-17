#include "Action_AvgBox.h"
#include "CpptrajStdio.h"
#include "DataSet_Mat3x3.h"

// Action_AvgBox::Help()
void Action_AvgBox::Help() const {
  mprintf("\t[<setname>] [out <file>]\n"
          "  Calculate average box.\n");
}

// Action_AvgBox::Init()
Action::RetType Action_AvgBox::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  std::string dsname = actionArgs.GetStringNext();

  // Add DataSet
  boxMatrix_ = init.DSL().AddSet( DataSet::MAT3X3, MetaData(dsname, "avg") );
  if (boxMatrix_ == 0) {
    mprinterr("Error: Could not allocate data set for average box.\n");
    return Action::ERR;
  }
# ifdef MPI
  // This set does not need to be synced since averaging is done by avgbox_
  boxMatrix_->SetNeedsSync( false );
# endif
  if (outfile != 0) outfile->AddDataSet( boxMatrix_ );

  mprintf("    AVGBOX: Calculating average box.\n");
  mprintf("\tAverage box set: %s\n", boxMatrix_->legend());
  if (outfile != 0) mprintf("\tOutput to file: %s\n", outfile->DataFilename().full());
  return Action::OK;
}

// Action_AvgBox::Setup()
Action::RetType Action_AvgBox::Setup(ActionSetup& setup)
{
  // Check box type
  if (!setup.CoordInfo().TrajBox().HasBox()) {
    mprintf("Error: Topology %s does not contain box information.\n",
            setup.Top().c_str());
    return Action::ERR;
  }
  return Action::OK;
}

// Action_AvgBox::DoAction()
Action::RetType Action_AvgBox::DoAction(int frameNum, ActionFrame& frm)
{
  Box const& box = frm.Frm().BoxCrd();
  Matrix_3x3 const& ucell = box.UnitCell();
  for (unsigned int i = 0; i < 9; i++)
    avgbox_[i].accumulate( ucell[i] );

  return Action::OK;
}

# ifdef MPI
int Action_AvgBox::SyncAction() {
  Matrix_3x3 ucell;
  // Transfer all processes avg to master and combine.
  if (trajComm_.Master()) {
    // Master
    double recvbuf[27];
    for (int proc = 1; proc < trajComm_.Size(); proc++) {
      trajComm_.Recv( recvbuf, 27, MPI_DOUBLE, proc, 2200 );
      int j = 0;
      for (int i = 0; i < 9; i++, j+= 3)
        avgbox_[i].Combine( Stats<double>(recvbuf[j], recvbuf[j+1], recvbuf[j+2]) );
    }
    for (int i = 0; i < 9; i++)
      ucell[i] = avgbox_[i].mean();
    ucell.Print("Average Unit Cell Vectors:");
  } else {
    // Children
    double sendbuf[27];
    int j = 0;
    for (int i = 0; i < 9; i++, j += 3) {
      sendbuf[j  ] = avgbox_[i].nData();
      sendbuf[j+1] = avgbox_[i].mean();
      sendbuf[j+2] = avgbox_[i].M2();
    }
    trajComm_.Send( sendbuf, 27, MPI_DOUBLE, 0, 2200 );
  }
  // Broadcast the unit cell
  ucell.BroadcastMatrix( trajComm_ );
  // Add to DataSet on all processes 
  ((DataSet_Mat3x3*)boxMatrix_)->AddMat3x3( ucell );
  return 0;
}
#endif

/** Do averaging in serial. In parallel just print box since averaging
  * is already done in SyncAction().
  */
void Action_AvgBox::Print() {
  mprintf("    AVGBOX:\n");
  if (avgbox_[0].nData() > 0) {
#   ifdef MPI
    // Unit cell data was already set up in SyncAction
    DataSet_Mat3x3 const& dset = static_cast<DataSet_Mat3x3 const&>( *boxMatrix_ );
    Matrix_3x3 const& ucell = dset[0];
#   else
    Matrix_3x3 ucell;
    for (int i = 0; i < 9; i++)
      ucell[i] = avgbox_[i].mean();
    ucell.Print("Average Unit Cell Vectors:");
    ((DataSet_Mat3x3*)boxMatrix_)->AddMat3x3( ucell );
#   endif
    Box thisBox;
    if (thisBox.SetupFromUcell( ucell )) {
      mprintf("Warning: Box appears to be invalid.\n");
    }
    thisBox.PrintInfo();
  } else {
    mprintf("\tNothing to average.\n");
  }
}
