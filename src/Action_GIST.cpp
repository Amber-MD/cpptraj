// GIST 
#include <cmath>
#include "Action_GIST.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG
#include "Action_Pairwise.h"

// CONSTRUCTOR
Action_GIST::Action_GIST() :
  gist_(0),
  watermodel_(false),
  useTIP3P_(false),
  useTIP4P_(false),
  useTIP4PEW_(false)
{
  gridcntr_[0] = -1;
  gridcntr_[1] = -1;
  gridcntr_[2] = -1;
  
  griddim_[0] = -1;
  griddim_[1] = -1;
  griddim_[2] = -1;
  
  gridspacn_ = 0;
 } 

void Action_GIST::Help() {
  mprintf("gist <watermodel>[{tip3p|tip4p|tip4pew}] [gridcntr <xval> <yval> <zval>] [griddim <xval> <yval> <zval>] [gridspacn <spaceval>] [out <filename>] \n");
  mprintf("\tCalculate GIST between water molecules in selected site \n");
}

// Action_GIST::init()
Action::RetType Action_GIST::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  //  InitImaging( !(actionArgs.hasKey("noimage")) );
  DataSet::scalarType stype = DataSet::UNDEFINED;
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  useTIP3P_ = actionArgs.hasKey("tip3p");
  useTIP4P_ = actionArgs.hasKey("tip4p");
  useTIP4PEW_ = actionArgs.hasKey("tip4pew");
  if (!useTIP3P_ && !useTIP4P_ && !useTIP4PEW_) {
    mprinterr("Error: gist: Only water models supprted are TIP3P and TIP4P\n");
    return Action::ERR;
  }

  if ( actionArgs.hasKey("gridcntr") ){
    gridcntr_[0] = actionArgs.getNextDouble(-1);
    gridcntr_[1] = actionArgs.getNextDouble(-1);
    gridcntr_[2] = actionArgs.getNextDouble(-1);
    mprintf("\tGIST grid center: %5.3f %5.3f %5.3f\n", gridcntr_[0],gridcntr_[1],gridcntr_[2]);
  }
  else{
    mprintf("\tGIST: No grid center values were found, using default\n");
    gridcntr_[0] = 0.0;
    gridcntr_[1] = 0.0;
    gridcntr_[2] = 0.0;
    mprintf("\tGIST grid center: %5.3f %5.3f %5.3f\n", gridcntr_[0],gridcntr_[1],gridcntr_[2]);
  }

  if ( actionArgs.hasKey("griddim") ){
    griddim_[0] = actionArgs.getNextDouble(-1);
    griddim_[1] = actionArgs.getNextDouble(-1);
    griddim_[2] = actionArgs.getNextDouble(-1);
    mprintf("\tGIST grid dimension: %5.3f %5.3f %5.3f\n", griddim_[0],griddim_[1],griddim_[2]);
  }
  else{
    mprintf("\tGIST: No grid dimensiom values were found, using default (box size) \n");
    //griddim_[0] = 0.0;
    //griddim_[1] = 0.0;
    //griddim_[2] = 0.0;
    //mprintf("\tGIST grid dimension: %5.3f %5.3f %5.3f\n", griddim_[0],griddim_[1],griddim_[2]);
  }

  gridspacn_ = actionArgs.getKeyDouble("gridspacn", 0.80);
  mprintf("\tGIST grid spacing: %5.3f \n", gridspacn_);

  // Dataset to store distances
  gist_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "Gist");
  if (gist_==0) return Action::ERR;
  gist_->SetScalar( DataSet::M_DISTANCE, stype );
  // Add dataset to data file
  if (outfile != 0) outfile->AddSet( gist_ );

  return Action::OK;
}

// Action_GIST::setup()
/** Set GIST up for this parmtop. Get masks etc.
  */
Action::RetType Action_GIST::Setup(Topology* currentParm, Topology** parmAddress) {
  mprintf("GIST Setup \n");

  return Action::OK;  
}


// Action_GIST::action()
Action::RetType Action_GIST::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {

  mprintf("GIST Action \n");
  //calculating energy
  return Action::OK;
}

