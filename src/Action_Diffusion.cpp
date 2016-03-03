#include <cmath> // sqrt
#include "Action_Diffusion.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // validDouble
#include "DataSet_1D.h" // LinearRegression
#ifdef TIMER
# include "Timer.h"
#endif

// CONSTRUCTOR
Action_Diffusion::Action_Diffusion() :
  avg_x_(0), avg_y_(0), avg_z_(0), avg_r_(0), avg_a_(0),
  printIndividual_(false),
  calcDiffConst_(false),
  time_(1.0),
  diffConst_(0),
  diffLabel_(0),
  diffSlope_(0),
  diffInter_(0),
  diffCorrl_(0),
  debug_(0),
  outputx_(0), outputy_(0), outputz_(0), outputr_(0), outputa_(0),
  diffout_(0),
  boxcenter_(0.0),
  masterDSL_(0)
{}

static inline void ShortHelp() {
  mprintf("\t[{out <filename> | separateout <suffix>}] [time <time per frame>]\n"
          "\t[<mask>] [<set name>] [individual] [diffout <filename>] [nocalc]\n");
}

void Action_Diffusion::Help() const {
  ShortHelp();
  mprintf("  Compute a mean square displacement plot for the atoms in the <mask>.\n"
          "  The following files are produced:\n"
          "    x_<suffix>: Mean square displacement(s) in the X direction (in Å^2).\n"
          "    y_<suffix>: Mean square displacement(s) in the Y direction (in Å^2).\n"
          "    z_<suffix>: Mean square displacement(s) in the Z direction (in Å^2).\n"
          "    r_<suffix>: Overall mean square displacement(s) (in Å^2).\n"
          "    a_<suffix>: Total distance travelled (in Å).\n");
}

static inline int CheckTimeArg(double dt) {
  if (dt <= 0.0) {
    mprinterr("Error: Diffusion time per frame incorrectly specified, must be > 0.0.\n");
    return 1;
  }
  return 0;
}

// Action_Diffusion::Init()
Action::RetType Action_Diffusion::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  debug_ = debugIn;
  image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  // Determine if this is old syntax or new.
  if (actionArgs.Nargs() > 2 && actionArgs.ArgIsMask(1) && validDouble(actionArgs[2]))
  {
    // Old syntax: <mask> <time per frame> [average] [<prefix>]
    printIndividual_ = !(actionArgs.hasKey("average"));
    calcDiffConst_ = false;
    mprintf("Warning: Deprecated syntax for 'diffusion'. Consider using new syntax:\n");
    ShortHelp();
    mask_.SetMaskString( actionArgs.GetMaskNext() );
    time_ = actionArgs.getNextDouble(1.0);
    if (CheckTimeArg(time_)) return Action::ERR;
    std::string outputNameRoot = actionArgs.GetStringNext();
    // Default filename: 'diffusion'
    if (outputNameRoot.empty())
      outputNameRoot.assign("diffusion");
    // Open output files
    ArgList oldArgs("prec 8.3 noheader");
    DataFile::DataFormatType dft = DataFile::DATAFILE;
    outputx_ = init.DFL().AddDataFile(outputNameRoot+"_x.xmgr", dft, oldArgs);
    outputy_ = init.DFL().AddDataFile(outputNameRoot+"_y.xmgr", dft, oldArgs);
    outputz_ = init.DFL().AddDataFile(outputNameRoot+"_z.xmgr", dft, oldArgs);
    outputr_ = init.DFL().AddDataFile(outputNameRoot+"_r.xmgr", dft, oldArgs);
    outputa_ = init.DFL().AddDataFile(outputNameRoot+"_a.xmgr", dft, oldArgs);
  } else {
    // New syntax: [{separateout <suffix> | out <filename>}] [time <time per frame>]
    //             [<mask>] [<set name>] [individual] [diffout <filename>] [nocalc]
    printIndividual_ = actionArgs.hasKey("individual");
    calcDiffConst_ = !(actionArgs.hasKey("nocalc"));
    std::string suffix = actionArgs.GetStringKey("separateout");
    std::string outname = actionArgs.GetStringKey("out");
    if (!outname.empty() && !suffix.empty()) {
      mprinterr("Error: Specify either 'out' or 'separateout', not both.\n");
      return Action::ERR;
    }
    diffout_ = init.DFL().AddDataFile(actionArgs.GetStringKey("diffout"));
    time_ = actionArgs.getKeyDouble("time", 1.0);
    if (CheckTimeArg(time_)) return Action::ERR;
    mask_.SetMaskString( actionArgs.GetMaskNext() );
    // Open output files.
    if (!suffix.empty()) {
      FileName FName( suffix );
      outputx_ = init.DFL().AddDataFile(FName.PrependFileName("x_"), actionArgs);
      outputy_ = init.DFL().AddDataFile(FName.PrependFileName("y_"), actionArgs);
      outputz_ = init.DFL().AddDataFile(FName.PrependFileName("z_"), actionArgs);
      outputr_ = init.DFL().AddDataFile(FName.PrependFileName("r_"), actionArgs);
      outputa_ = init.DFL().AddDataFile(FName.PrependFileName("a_"), actionArgs);
      if (diffout_ == 0 && calcDiffConst_)
        diffout_ = init.DFL().AddDataFile(FName.PrependFileName("diff_"), actionArgs);
    } else if (!outname.empty()) {
      outputr_ = init.DFL().AddDataFile( outname, actionArgs );
      outputx_ = outputy_ = outputz_ = outputa_ = outputr_;
    }
  }
  if (diffout_ != 0) calcDiffConst_ = true;
  // Add DataSets
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
  if (outputr_ != 0) outputr_->AddDataSet( avg_r_ );
  if (outputx_ != 0) outputx_->AddDataSet( avg_x_ );
  if (outputy_ != 0) outputy_->AddDataSet( avg_y_ );
  if (outputz_ != 0) outputz_->AddDataSet( avg_z_ );
  if (outputa_ != 0) outputa_->AddDataSet( avg_a_ );
  // Set X dim
  Xdim_ = Dimension(0.0, time_, "Time");
  avg_x_->SetDim(Dimension::X, Xdim_);
  avg_y_->SetDim(Dimension::X, Xdim_);
  avg_z_->SetDim(Dimension::X, Xdim_);
  avg_r_->SetDim(Dimension::X, Xdim_);
  avg_a_->SetDim(Dimension::X, Xdim_);
  // Add DataSets for diffusion constant calc
  if (calcDiffConst_) {
    MetaData::tsType ts = MetaData::NOT_TS;
    diffConst_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "D", ts));
    diffLabel_ = init.DSL().AddSet(DataSet::STRING, MetaData(dsname_, "Label", ts));
    diffSlope_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "Slope", ts));
    diffInter_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "Intercept", ts));
    diffCorrl_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "Corr", ts));
    if (diffConst_ == 0 || diffLabel_ == 0 || diffSlope_ == 0 || diffInter_ == 0 ||
        diffCorrl_ == 0)
      return Action::ERR;
#   ifdef MPI
    // No sync needed since these are not time series
    diffConst_->SetNeedsSync( false );
    diffLabel_->SetNeedsSync( false );
    diffSlope_->SetNeedsSync( false );
    diffInter_->SetNeedsSync( false );
    diffCorrl_->SetNeedsSync( false );
#   endif
    if (diffout_ != 0) {
      diffout_->AddDataSet( diffConst_ );
      diffout_->AddDataSet( diffSlope_ );
      diffout_->AddDataSet( diffInter_ );
      diffout_->AddDataSet( diffCorrl_ );
      diffout_->AddDataSet( diffLabel_ );
    }
    Dimension Ddim( 1, 1, "Set" );
    diffConst_->SetDim(Dimension::X, Ddim);
    diffLabel_->SetDim(Dimension::X, Ddim);
    diffSlope_->SetDim(Dimension::X, Ddim);
    diffInter_->SetDim(Dimension::X, Ddim);
    diffCorrl_->SetDim(Dimension::X, Ddim);
  }
  // Save master data set list, needed when printIndividual_
  masterDSL_ = init.DslPtr();

  mprintf("    DIFFUSION:\n");
  mprintf("\tAtom Mask is [%s]\n", mask_.MaskString());
  if (printIndividual_)
    mprintf("\tBoth average and individual diffusion will be calculated.\n");
  else
    mprintf("\tOnly average diffusion will be calculated.\n");
  mprintf("\tData set base name: %s\n", avg_x_->Meta().Name().c_str());
  if (image_.UseImage())
    mprintf("\tCorrections for imaging enabled.\n");
  else
    mprintf("\tCorrections for imaging disabled.\n");
  // If one file defined, assume all are.
  if (outputx_ != 0) {
    mprintf("\tOutput files:\n"
            "\t  %s: (x) Mean square displacement(s) in the X direction (in Å^2).\n"
            "\t  %s: (y) Mean square displacement(s) in the Y direction (in Å^2).\n"
            "\t  %s: (z) Mean square displacement(s) in the Z direction (in Å^2).\n"
            "\t  %s: (r) Overall mean square displacement(s) (in Å^2).\n"
            "\t  %s: (a) Total distance travelled (in Å).\n",
            outputx_->DataFilename().full(), outputy_->DataFilename().full(),
            outputz_->DataFilename().full(), outputr_->DataFilename().full(),
            outputa_->DataFilename().full());
  }
  mprintf("\tThe time between frames is %g ps.\n", time_);
  if (calcDiffConst_) {
    mprintf("\tCalculating diffusion constants by fitting slope to MSD vs time\n"
            "\t  and multiplying by 10.0/2*N (where N is # of dimensions), units\n"
            "\t  are 1x10^-5 cm^2/s.\n");
    if (diffout_ != 0)
      mprintf("\tDiffusion constant output to '%s'\n", diffout_->DataFilename().full());
    else
      mprintf("\tDiffusion constant output to STDOUT.\n");
  } else
    mprintf("\tTo calculate diffusion constant from mean squared displacement plots,\n"
            "\t  calculate the slope of MSD vs time and multiply by 10.0/2*N (where N\n"
            "\t  is # of dimensions); this will give units of 1x10^-5 cm^2/s.\n");
  return Action::OK;
}

// Action_Diffusion::Setup()
Action::RetType Action_Diffusion::Setup(ActionSetup& setup) {
# ifdef TIMER
  Timer time_setup, time_addsets;
  time_setup.Start();
# endif
  // Setup atom mask
  if (setup.Top().SetupIntegerMask( mask_ )) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: No atoms selected.\n");
    return Action::SKIP;
  }

  // Set up imaging info for this parm
  image_.SetupImaging( setup.CoordInfo().TrajBox().Type() );
  if (image_.ImagingEnabled())
    mprintf("\tImaging enabled.\n");
  else
    mprintf("\tImaging disabled.\n");

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
    // Create as many spots for sets as needed. All do not have to be used.
    if (mask_.back() >= (int)atom_x_.size()) {
      int newSize = mask_.back() + 1;
      atom_x_.resize( newSize, 0 );
      atom_y_.resize( newSize, 0 );
      atom_z_.resize( newSize, 0 );
      atom_r_.resize( newSize, 0 );
      atom_a_.resize( newSize, 0 );
    }
    for (AtomMask::const_iterator at = mask_.begin(); at != mask_.end(); at++)
    {
      if (atom_x_[*at] == 0) {
#       ifdef TIMER
        time_addsets.Start();
#       endif
        atom_x_[*at] = masterDSL_->AddSet_NoCheck(DataSet::FLOAT, MetaData(dsname_, "aX", *at+1));
        atom_y_[*at] = masterDSL_->AddSet_NoCheck(DataSet::FLOAT, MetaData(dsname_, "aY", *at+1));
        atom_z_[*at] = masterDSL_->AddSet_NoCheck(DataSet::FLOAT, MetaData(dsname_, "aZ", *at+1));
        atom_r_[*at] = masterDSL_->AddSet_NoCheck(DataSet::FLOAT, MetaData(dsname_, "aR", *at+1));
        atom_a_[*at] = masterDSL_->AddSet_NoCheck(DataSet::FLOAT, MetaData(dsname_, "aA", *at+1));
#       ifdef TIMER
        time_addsets.Stop();
#       endif
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
# ifdef TIMER
  time_setup.Stop();
  time_addsets.WriteTiming(3, "Diffusion Add Sets", time_setup.Total());
  time_setup.WriteTiming(2, "Diffusion Setup");
# endif
  return Action::OK;
}

// Action_Diffusion::DoAction()
Action::RetType Action_Diffusion::DoAction(int frameNum, ActionFrame& frm) {
  Matrix_3x3 ucell, recip;
  // Load initial frame if necessary
  if (initial_.empty()) {
    initial_ = frm.Frm();
#   ifdef MPI
    if (trajComm_.Size() > 1) {
      if (trajComm_.Master())
        for (int rank = 1; rank < trajComm_.Size(); ++rank)
          initial_.SendFrame( rank, trajComm_ );
      else
        initial_.RecvFrame( 0, trajComm_ );
    }
#   endif
    for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
    {
      const double* XYZ = initial_.XYZ(*atom); //frm.Frm().XYZ(*atom);
      previous_.push_back( XYZ[0] );
      previous_.push_back( XYZ[1] );
      previous_.push_back( XYZ[2] );
    }
#   ifdef MPI
    if (trajComm_.Master()) return Action::OK;
    // On non-master we want to calculate diffusion for *this* frame. It is
    // possible wrapping has already occurred from the initial frame. Try to
    // detect and correct for this.
#   endif
  }
  // Diffusion calculation
  if (image_.ImageType() != NOIMAGE) {
    boxcenter_ = frm.Frm().BoxCrd().Center();
    if (image_.ImageType() == NONORTHO)
      frm.Frm().BoxCrd().ToRecip(ucell, recip);
  }
  // For averaging over selected atoms
  double average2 = 0.0;
  double avgx = 0.0;
  double avgy = 0.0;
  double avgz = 0.0;
  unsigned int idx = 0; // Index into previous_ and delta_
  mprintf("\nDEBUG: Diffusion (%i): boxcenter={ %g %g %g }\n", frameNum+1, boxcenter_[0], boxcenter_[1], boxcenter_[2]);
  for (AtomMask::const_iterator at = mask_.begin(); at != mask_.end(); ++at, idx += 3)
  { // Get current and initial coords for this atom.
    const double* XYZ = frm.Frm().XYZ(*at);
    const double* iXYZ = initial_.XYZ(*at);
    mprintf("\tInitial={ %g %g %g }  Current={ %g %g %g }\n",
            iXYZ[0], iXYZ[1], iXYZ[2], XYZ[0], XYZ[1], XYZ[2]);
    // Calculate distance to previous frames coordinates.
    double delx = XYZ[0] - previous_[idx  ];
    double dely = XYZ[1] - previous_[idx+1];
    double delz = XYZ[2] - previous_[idx+2];
    mprintf("\tCurrentDeltaFromPrevious={ %g %g %g }\n", delx, dely, delz); // DEBUG
    if (image_.ImageType() == NOIMAGE) {
      // No imaging. Calculate distance from current position to initial position.
      delx = XYZ[0] - iXYZ[0];
      dely = XYZ[1] - iXYZ[1];
      delz = XYZ[2] - iXYZ[2];
    } else if ( image_.ImageType() == ORTHO ) {
      // If the particle moved more than half the box, assume it was imaged
      // and adjust the distance of the total movement with respect to the
      // original frame.
      if      (delx >  boxcenter_[0]) delta_[idx  ] -= frm.Frm().BoxCrd().BoxX();
      else if (delx < -boxcenter_[0]) delta_[idx  ] += frm.Frm().BoxCrd().BoxX();
      if      (dely >  boxcenter_[1]) delta_[idx+1] -= frm.Frm().BoxCrd().BoxY();
      else if (dely < -boxcenter_[1]) delta_[idx+1] += frm.Frm().BoxCrd().BoxY();
      if      (delz >  boxcenter_[2]) delta_[idx+2] -= frm.Frm().BoxCrd().BoxZ();
      else if (delz < -boxcenter_[2]) delta_[idx+2] += frm.Frm().BoxCrd().BoxZ();
      //mprinterr("\tDelta={ %g %g %g }\n", delta_[idx], delta_[idx+1], delta_[idx+2]); // DEBUG
      // Calculate the distance between this "fixed" coordinate
      // and the reference (initial) frame.
      delx = XYZ[0] + delta_[idx  ] - iXYZ[0];
      dely = XYZ[1] + delta_[idx+1] - iXYZ[1];
      delz = XYZ[2] + delta_[idx+2] - iXYZ[2];
    } else if ( image_.ImageType() == NONORTHO ) {
      // Non-orthorhombic imaging
      if (fabs(delx) > boxcenter_[0] ||
          fabs(dely) > boxcenter_[1] ||
          fabs(delz) > boxcenter_[2])
      {
        mprintf("\tImaging detected.\n");
        // Previous position in fractional coords
        //Vec3 pFrac = recip * Vec3( previous_[idx], previous_[idx+1], previous_[idx+2] );
        // Previous position back in Cartesian space
        //Vec3 pCart = ucell.TransposeMult( pFrac );
        Vec3 pCart( previous_[idx], previous_[idx+1], previous_[idx+2] );
        //mprintf("\tPrevious position (Cart): %g %g %g (%g %g %g)\n",
        //        previous_[idx], previous_[idx+1], previous_[idx+2], pCart[0], pCart[1], pCart[2]);
        // Current position in fractional coords
        Vec3 cFrac = recip * Vec3( XYZ[0], XYZ[1], XYZ[2] );
        // Current position back in Cartesian space
        //Vec3 cCart = ucell.TransposeMult( cFrac );
        // Look for imaged distance closer than current position
        double minDist2 = frm.Frm().BoxCrd().BoxX() *
                          frm.Frm().BoxCrd().BoxY() *
                          frm.Frm().BoxCrd().BoxZ();
        Vec3 minCurr(0.0);
        for (int ix = -1; ix < 2; ix++) {
          for (int iy = -1; iy < 2; iy++) {
            for (int iz = -1; iz < 2; iz++) {
              if (ix != 0 || iy != 0 || iz != 0) { // Ignore current position
                Vec3 ixyz(ix, iy, iz);
                // Current position shifted and back in Cartesian space
                Vec3 IMG = ucell.TransposeMult(cFrac + ixyz);
                // Distance from previous position to imaged current position
                Vec3 dxyz = IMG - pCart;
                double dist2 = dxyz.Magnitude2();
                if (dist2 < minDist2) {
                  minDist2 = dist2;
                  minCurr = IMG;
                }
              }
            }
          }
        }
        mprintf("\tClosest position {%g %g %g} (%g) to previous={ %g %g %g }\n",
                minCurr[0], minCurr[1], minCurr[2], sqrt(minDist2), pCart[0], pCart[1], pCart[2]);
                
        // Update the delta for this atom
        delta_[idx  ] += minCurr[0] - XYZ[0]; // cCart
        delta_[idx+1] += minCurr[1] - XYZ[1];
        delta_[idx+2] += minCurr[2] - XYZ[2];
        mprintf("\tUpdated delta={ %g %g %g }\n", delta_[idx], delta_[idx+1], delta_[idx+2]);
      }
      // Calculate the distance between this "fixed" coordinate
      // and the reference (initial) frame.
      delx = XYZ[0] + delta_[idx  ] - iXYZ[0];
      dely = XYZ[1] + delta_[idx+1] - iXYZ[1];
      delz = XYZ[2] + delta_[idx+2] - iXYZ[2];
      mprintf("\tCorrected distance from initial={ %g %g %g }\n", delx, dely, delz);
    }
    //mprinterr("\tDeltaFromInitial={ %g %g %g }\n", delx, dely, delz); // DEBUG
    // Calc distances for this atom
    double distx = delx * delx;
    double disty = dely * dely;
    double distz = delz * delz;
    double dist2 = distx + disty + distz;
    mprintf("\tDistances(X Y Z Tot): %g %g %g %g\n", sqrt(distx), sqrt(disty), sqrt(distz), sqrt(dist2));
    // Accumulate averages
    avgx += distx;
    avgy += disty;
    avgz += distz;
    average2 += dist2;
    // Store distances for this atom
    if (printIndividual_) {
      float fval = (float)distx;
      atom_x_[*at]->Add(frameNum, &fval);
      fval = (float)disty;
      atom_y_[*at]->Add(frameNum, &fval);
      fval = (float)distz;
      atom_z_[*at]->Add(frameNum, &fval);
      fval = (float)dist2;
      atom_r_[*at]->Add(frameNum, &fval);
      dist2 = sqrt(dist2);
      fval = (float)dist2;
      atom_a_[*at]->Add(frameNum, &fval);
    }
    // Update the previous coordinate set to match the current coordinates
    previous_[idx  ] = XYZ[0];
    previous_[idx+1] = XYZ[1];
    previous_[idx+2] = XYZ[2];
  } // END loop over selected atoms
  // Calc averages
  double dNselected = 1.0 / (double)mask_.Nselected();
  avgx *= dNselected;
  avgy *= dNselected;
  avgz *= dNselected;
  average2 *= dNselected;
  // Save averages
  avg_x_->Add(frameNum, &avgx);
  avg_y_->Add(frameNum, &avgy);
  avg_z_->Add(frameNum, &avgz);
  avg_r_->Add(frameNum, &average2);
  average2 = sqrt(average2);
  avg_a_->Add(frameNum, &average2);
  return Action::OK;
}

void Action_Diffusion::Print() {
  if (!calcDiffConst_) return;
  mprintf("    DIFFUSION: Calculating diffusion constants from slopes.\n");
  std::string const& name = avg_r_->Meta().Name();
  unsigned int set = 0;
  CalcDiffusionConst( set, avg_r_, 3, name + "_AvgDr" );
  CalcDiffusionConst( set, avg_x_, 1, name + "_AvgDx" );
  CalcDiffusionConst( set, avg_y_, 1, name + "_AvgDy" );
  CalcDiffusionConst( set, avg_z_, 1, name + "_AvgDz" );
  if (printIndividual_) {
    CalcDiffForSet( set, atom_r_, 3, name + "_dr" );
    CalcDiffForSet( set, atom_x_, 3, name + "_dx" );
    CalcDiffForSet( set, atom_y_, 3, name + "_dy" );
    CalcDiffForSet( set, atom_z_, 3, name + "_dz" );
  }
}

void Action_Diffusion::CalcDiffForSet(unsigned int& set, Dlist const& Sets, int Ndim,
                                      std::string const& label) const
{
  for (Dlist::const_iterator ds = Sets.begin(); ds != Sets.end(); ds++)
    if (*ds != 0)
      CalcDiffusionConst(set, *ds, Ndim, label + "_" + integerToString( (*ds)->Meta().Idx() ));
}

void Action_Diffusion::CalcDiffusionConst(unsigned int& set, DataSet* ds, int Ndim,
                                          std::string const& label) const
{
  DataSet_1D const& data = static_cast<DataSet_1D const&>( *ds );
  double Factor = 10.0 / ((double)Ndim * 2.0);
  double slope, intercept, corr;
  double Dval = 0.0;
  if (data.LinearRegression( slope, intercept, corr, 0 ) == 0)
    Dval = slope * Factor;
  if (diffout_ == 0)
    mprintf("\t'%s' D= %g  Slope= %g  Int= %g  Corr= %g\n", data.legend(), Dval,
            slope, intercept, corr);
  diffConst_->Add(set  , &Dval);
  diffSlope_->Add(set  , &slope);
  diffInter_->Add(set  , &intercept);
  diffCorrl_->Add(set  , &corr);
  diffLabel_->Add(set++, label.c_str());
}
