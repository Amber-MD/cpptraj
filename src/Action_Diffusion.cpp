#include <cmath> // sqrt
#include "Action_Diffusion.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Diffusion::Action_Diffusion() :
  avg_x_(0), avg_y_(0), avg_z_(0), avg_r_(0), avg_a_(0),
  printIndividual_(false),
  time_(1.0),
  hasBox_(false),
  debug_(0),
  outputx_(0), outputy_(0), outputz_(0), outputr_(0), outputa_(0),
  boxcenter_(0.0),
  masterDSL_(0)
{}

void Action_Diffusion::Help() {
  mprintf("\t<mask> <time per frame> [average] [<prefix>]\n"
          "  Compute a mean square displacement plot for the atoms in the mask.\n"
          "  The following files are produced:\n"
          "    <prefix>_x.xmgr: Mean square displacement(s) in the X direction (in Å^2).\n"
          "    <prefix>_y.xmgr: Mean square displacement(s) in the Y direction (in Å^2).\n"
          "    <prefix>_z.xmgr: Mean square displacement(s) in the Z direction (in Å^2).\n"
          "    <prefix>_r.xmgr: Overall mean square displacement(s) (in Å^2).\n"
          "    <prefix>_a.xmgr: Total distance travelled (in Å).\n");
}

Action::RetType Action_Diffusion::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  printIndividual_ = !(actionArgs.hasKey("average"));
  mask_.SetMaskString( actionArgs.GetMaskNext() );
  time_ = actionArgs.getNextDouble(1.0);
  if (time_ < 0) {
    mprinterr("Error: diffusion: time per frame incorrectly specified\n");
    return Action::ERR;
  }
  std::string outputNameRoot = actionArgs.GetStringNext();
  // Default filename: 'diffusion'
  if (outputNameRoot.empty()) 
    outputNameRoot.assign("diffusion");
  
  // Open output files
  outputx_ = init.DFL().AddDataFile(outputNameRoot+"_x.xmgr");//, "X MSD");
  outputy_ = init.DFL().AddDataFile(outputNameRoot+"_y.xmgr");//, "Y MSD");
  outputz_ = init.DFL().AddDataFile(outputNameRoot+"_z.xmgr");//, "Z MSD");
  outputr_ = init.DFL().AddDataFile(outputNameRoot+"_r.xmgr");//, "Overall MSD");
  outputa_ = init.DFL().AddDataFile(outputNameRoot+"_a.xmgr");//, "Total Distance");
  //if (outputx_ == 0 || outputy_ == 0 || outputz_ == 0 ||
  //    outputr_ == 0 || outputa_ == 0) return Action::ERR;
  // AddData Sets
  dsname_ = actionArgs.GetStringNext();
  if (dsname_.empty())
    dsname_ = init.DSL().GenerateDefaultName("Diff");
  avg_x_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "X"));
  avg_y_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "Y"));
  avg_z_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "Z"));
  avg_r_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "R"));
  avg_a_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "A"));
  if (avg_x_ == 0 || avg_y_ == 0 || avg_z_ == 0 || avg_r_ == 0 || avg_a_ == 0)
    return Action::ERR;
  if (outputx_ != 0) outputx_->AddDataSet( avg_x_ );
  if (outputy_ != 0) outputy_->AddDataSet( avg_y_ );
  if (outputz_ != 0) outputz_->AddDataSet( avg_z_ );
  if (outputr_ != 0) outputr_->AddDataSet( avg_r_ );
  if (outputa_ != 0) outputa_->AddDataSet( avg_a_ );
  // Set X dim
  Xdim_ = Dimension(0.0, time_);
  avg_x_->SetDim(Dimension::X, Xdim_);
  avg_y_->SetDim(Dimension::X, Xdim_);
  avg_z_->SetDim(Dimension::X, Xdim_);
  avg_r_->SetDim(Dimension::X, Xdim_);
  avg_a_->SetDim(Dimension::X, Xdim_);
  // Save master data set list
  masterDSL_ = init.DslPtr();

  mprintf("    DIFFUSION:\n");
  mprintf("\tAtom Mask is [%s]\n", mask_.MaskString());
  if (printIndividual_)
    mprintf("\tThe average and individual results will be printed to:\n");
  else
    mprintf("\tOnly the average results will be printed to:\n");
  const char* onr = outputNameRoot.c_str();
  mprintf("\t  %s_x.xmgr: Mean square displacement(s) in the X direction (in Å^2).\n"
          "\t  %s_y.xmgr: Mean square displacement(s) in the Y direction (in Å^2).\n"
          "\t  %s_z.xmgr: Mean square displacement(s) in the Z direction (in Å^2).\n"
          "\t  %s_r.xmgr: Overall mean square displacement(s) (in Å^2).\n"
          "\t  %s_a.xmgr: Total distance travelled (in Å).\n",
          onr, onr, onr, onr, onr);
  mprintf("\tThe time between frames in ps is %.3f.\n", time_);
  mprintf("\tTo calculate diffusion constants from a mean squared displacement plot\n"
          "\t(i.e. {_x|_y|_z|_r}.xmgr), calculate the slope of the line and multiply\n"
          "\tby 10.0/6.0; this will give units of 1x10^-5 cm^2/s\n");

  return Action::OK;
}

Action::RetType Action_Diffusion::Setup(ActionSetup& setup) {
  // Setup atom mask
  if (setup.Top().SetupIntegerMask( mask_ )) return Action::ERR;
  if (mask_.None()) {
    mprintf("Warning: No atoms selected.\n");
    return Action::SKIP;
  }

  // Check for box
  if ( setup.CoordInfo().TrajBox().Type() != Box::NOBOX ) {
    // Currently only works for orthogonal boxes
    if ( setup.CoordInfo().TrajBox().Type() != Box::ORTHO ) {
      mprintf("Warning: diffusion command currently only works with orthogonal boxes.\n");
      mprintf("Warning: imaging will be disabled for this command. This may result in\n");
      mprintf("Warning: large jumps if target molecule is imaged. To counter this please\n");
      mprintf("Warning: use the 'unwrap' command prior to 'diffusion'.\n");
      hasBox_ = false;
    } else 
      hasBox_ = true;
  } else
    hasBox_ = false;

  // Allocate the delta array
  delta_.assign( mask_.Nselected() * 3, 0.0 );

  // Reserve space for the previous coordinates array
  previous_.reserve( mask_.Nselected() * 3 );

  // If initial frame already set and current # atoms > # atoms in initial
  // frame this will probably cause an error.
  if (!initial_.empty() && setup.Top().Natom() > initial_.Natom()) {
    mprintf("Warning: # atoms in current parm (%s, %i) > # atoms in initial frame (%i)\n",
             setup.Top().c_str(), setup.Top().Natom(), initial_.Natom());
    mprintf("Warning: This may lead to segmentation faults.\n");
  }

  // Set up sets for individual atoms if necessary
  if (printIndividual_) {
    for (AtomMask::const_iterator at = mask_.begin(); at != mask_.end(); at++)
    {
      if (*at >= (int)atom_x_.size()) {
        int newSize = *at + 1;
        atom_x_.resize( newSize, 0 );
        atom_y_.resize( newSize, 0 );
        atom_z_.resize( newSize, 0 );
        atom_r_.resize( newSize, 0 );
        atom_a_.resize( newSize, 0 );
      }
      if (atom_x_[*at] == 0) { // TODO: FLOAT?
        atom_x_[*at] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "aX", *at+1));
        atom_y_[*at] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "aY", *at+1));
        atom_z_[*at] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "aZ", *at+1));
        atom_r_[*at] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "aR", *at+1));
        atom_a_[*at] = masterDSL_->AddSet(DataSet::DOUBLE, MetaData(dsname_, "aA", *at+1));
        if (outputx_ != 0) outputx_->AddDataSet(atom_x_[*at]);
        if (outputy_ != 0) outputy_->AddDataSet(atom_y_[*at]);
        if (outputz_ != 0) outputz_->AddDataSet(atom_z_[*at]);
        if (outputr_ != 0) outputr_->AddDataSet(atom_r_[*at]);
        if (outputa_ != 0) outputa_->AddDataSet(atom_a_[*at]);
        atom_x_[*at]->SetDim(Dimension::X, Xdim_);
        atom_y_[*at]->SetDim(Dimension::X, Xdim_);
        atom_z_[*at]->SetDim(Dimension::X, Xdim_);
        atom_r_[*at]->SetDim(Dimension::X, Xdim_);
        atom_a_[*at]->SetDim(Dimension::X, Xdim_);
      }
    }
  }

  return Action::OK;
}

Action::RetType Action_Diffusion::DoAction(int frameNum, ActionFrame& frm) {
  // Load initial frame if necessary
  if (initial_.empty()) {
    initial_ = frm.Frm();
    for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
    {
      const double* XYZ = frm.Frm().XYZ(*atom);
      previous_.push_back( XYZ[0] );
      previous_.push_back( XYZ[1] );
      previous_.push_back( XYZ[2] );
    }
  } else {
    if (hasBox_) 
      boxcenter_ = frm.Frm().BoxCrd().Center();
    Vec3 boxL = frm.Frm().BoxCrd().Lengths();
    // For averaging over selected atoms
    double average2 = 0.0;
    double avgx = 0.0;
    double avgy = 0.0;
    double avgz = 0.0;
    unsigned int idx = 0;
    for (AtomMask::const_iterator at = mask_.begin(); at != mask_.end(); ++at, idx += 3)
    { // Get current and initial coords for this atom.
      const double* XYZ = frm.Frm().XYZ(*at);
      const double* iXYZ = initial_.XYZ(*at);
      // Calculate distance to previous frames coordinates.
      double delx = XYZ[0] - previous_[idx  ];
      double dely = XYZ[1] - previous_[idx+1];
      double delz = XYZ[2] - previous_[idx+2];
      // If the particle moved more than half the box, assume
      // it was imaged and adjust the distance of the total
      // movement with respect to the original frame.
      if (hasBox_) {
        if      (delx >  boxcenter_[0]) delta_[idx  ] -= boxL[0];
        else if (delx < -boxcenter_[0]) delta_[idx  ] += boxL[0];
        else if (dely >  boxcenter_[1]) delta_[idx+1] -= boxL[1];
        else if (dely < -boxcenter_[1]) delta_[idx+1] += boxL[1];
        else if (delz >  boxcenter_[2]) delta_[idx+2] -= boxL[2];
        else if (delz < -boxcenter_[2]) delta_[idx+2] += boxL[2];
      }
      // DEBUG
      if (debug_ > 2)
        mprintf("ATOM: %5i %10.3f %10.3f %10.3f",*at,XYZ[0],delx,delta_[idx  ]);
      // Set the current x with reference to the un-imaged trajectory.
      double xx = XYZ[0] + delta_[idx  ]; 
      double yy = XYZ[1] + delta_[idx+1]; 
      double zz = XYZ[2] + delta_[idx+2];
      // Calculate the distance between this "fixed" coordinate
      // and the reference (initial) frame.
      delx = xx - iXYZ[0];
      dely = yy - iXYZ[1];
      delz = zz - iXYZ[2];
      // DEBUG
      if (debug_ > 2)
        mprintf(" %10.3f\n", delx);
      // Calc distances for this atom
      double distx = delx * delx;
      double disty = dely * dely;
      double distz = delz * delz;
      double dist2 = distx + disty + distz;
      // Accumulate averages
      avgx += distx;
      avgy += disty;
      avgz += distz;
      average2 += dist2;
      // Store distances for this atom
      if (printIndividual_) {
        atom_x_[*at]->Add(frameNum, &distx);
        atom_y_[*at]->Add(frameNum, &disty);
        atom_z_[*at]->Add(frameNum, &distz);
        atom_r_[*at]->Add(frameNum, &dist2);
        dist2 = sqrt(dist2);
        atom_a_[*at]->Add(frameNum, &dist2);
      }
      // Update the previous coordinate set to match the current coordinates
      previous_[idx  ] = XYZ[0];
      previous_[idx+1] = XYZ[1];
      previous_[idx+2] = XYZ[2];
    } // END loop over selected atoms
    double dNselected = (double)mask_.Nselected();
    avgx /= dNselected;
    avgy /= dNselected;
    avgz /= dNselected;
    average2 /= dNselected;

    avg_x_->Add(frameNum, &avgx);
    avg_y_->Add(frameNum, &avgy);
    avg_z_->Add(frameNum, &avgz);
    avg_r_->Add(frameNum, &average2);
    average2 = sqrt(average2);
    avg_a_->Add(frameNum, &average2);
/*
    // ----- OUTPUT -----
    // Output averages
    double Time = time_ * (double)frameNum;
    outputx_->Printf("%8.3f  %8.3f", Time, avgx);
    outputy_->Printf("%8.3f  %8.3f", Time, avgy);
    outputz_->Printf("%8.3f  %8.3f", Time, avgz);
    outputr_->Printf("%8.3f  %8.3f", Time, average);
    outputa_->Printf("%8.3f  %8.3f", Time, sqrt(average));
    // Individual values
    if (printIndividual_) {
      for (int i = 0; i < mask_.Nselected(); ++i) {
        outputx_->Printf("  %8.3f", distancex_[i]);
        outputy_->Printf("  %8.3f", distancey_[i]);
        outputz_->Printf("  %8.3f", distancez_[i]);
        outputr_->Printf("  %8.3f", distance_[i]);
        outputa_->Printf("  %8.3f", sqrt(distance_[i]));
      }
    }
    // Print newlines
    outputx_->Printf("\n");
    outputy_->Printf("\n");
    outputz_->Printf("\n");
    outputr_->Printf("\n");
    outputa_->Printf("\n");
*/
  }
  return Action::OK;
}  
