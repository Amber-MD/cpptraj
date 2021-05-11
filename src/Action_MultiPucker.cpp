#include "Action_MultiPucker.h"
#include "CpptrajStdio.h"
#include "TorsionRoutines.h"

using namespace Cpptraj;

const double Action_MultiPucker::PERIOD_ = 360.0;

/** CONSTRUCTOR */
Action_MultiPucker::Action_MultiPucker() :
  outfile_(0),
  ampfile_(0),
  thetafile_(0),
  masterDSL_(0),
  defaultMethod_(Pucker::ALTONA_SUNDARALINGAM),
  puckerMin_(0),
  puckerMax_(0),
  offset_(0),
  calc_amp_(false),
  calc_theta_(false)
{}

// Action_MultiPucker::Help()
void Action_MultiPucker::Help() const {
  mprintf("\t[<name>] [<pucker types>] [resrange <range>] [out <filename>]\n"
          "\t[altona|cremer] [%s]\n"
          "\t[amplitude [ampout <ampfile>]] [theta [thetaout <thetafile>]]\n"
          "\t[range360] [offset <offset>]\n",
          Pucker::PuckerSearch::newTypeArgsHelp());
  mprintf("\t<pucker types> = ");
  Pucker::PuckerSearch::ListKnownTypes();
  mprintf("  Calculate specified pucker types for residues in given <range>.\n");
}

// Action_MultiPucker::Init()
Action::RetType Action_MultiPucker::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get Keywords
  outfile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs);
  if      (actionArgs.hasKey("altona")) defaultMethod_ = Pucker::ALTONA_SUNDARALINGAM;
  else if (actionArgs.hasKey("cremer")) defaultMethod_ = Pucker::CREMER_POPLE;
  else                                  defaultMethod_ = Pucker::UNSPECIFIED;
  ampfile_ = 0;
  calc_amp_ = actionArgs.hasKey("amplitude");
  if (calc_amp_) {
    ampfile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("ampout"), actionArgs);
  }
  calc_theta_ = actionArgs.hasKey("theta");
  if (calc_theta_) {
    thetafile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("thetaout"), actionArgs);
  }
  offset_ = actionArgs.getKeyDouble("offset",0.0);
  if (actionArgs.hasKey("range360"))
    puckerMin_ = 0.0;
  else
    puckerMin_ = -180.0;
  puckerMax_ = puckerMin_ + PERIOD_;
  std::string resrange_arg = actionArgs.GetStringKey("resrange");
  if (!resrange_arg.empty())
    if (resRange_.SetRange( resrange_arg )) return Action::ERR;
  // Search for known pucker keywords
  if (puckerSearch_.SearchForArgs( actionArgs )) return Action::ERR;
  // Get custom pucker args
  if (puckerSearch_.SearchForNewTypeArgs( actionArgs )) return Action::ERR;
  // If no pucker types are yet selected, this will select all.
  puckerSearch_.SearchForAll();

  // Setup DataSet(s) name
  dsetname_ = actionArgs.GetStringNext();

  mprintf("    MULTIPUCKER: Calculating");
  puckerSearch_.PrintTypes();
  if (!resRange_.Empty())
    mprintf(" puckers for residues in range %s\n", resRange_.RangeArg());
  else
    mprintf(" puckers for all solute residues.\n");
  if (defaultMethod_ == Pucker::ALTONA_SUNDARALINGAM) 
    mprintf("\tUsing Altona & Sundaralingam method.\n");
  else if (defaultMethod_ == Pucker::CREMER_POPLE)
    mprintf("\tUsing Cremer & Pople method.\n");
  else
    mprintf("\tAltona & Sundaralingam method will be used for 5 atom puckers.\n"
            "\tCremer & Pople method will be used for 6 atom puckers.\n");
  if (offset_!=0)
    mprintf("\tOffset: %f degrees will be added to values.\n", offset_);
  if (puckerMin_ > -180.0)
    mprintf("\tOutput range is 0 to 360 degrees.\n");
  else
    mprintf("\tOutput range is -180 to 180 degrees.\n");
  if (!dsetname_.empty())
    mprintf("\tDataSet name: %s\n", dsetname_.c_str());
  if (outfile_ != 0) mprintf("\tOutput to %s\n", outfile_->DataFilename().full());
  if (calc_amp_) {
    mprintf("\tAmplitudes (in degrees) will be calculated.\n");
    if (ampfile_ != 0) mprintf("\tAmplitudes output to %s\n", ampfile_->DataFilename().full());
  }
  if (calc_theta_) {
    mprintf("\tThetas (in degrees) will be calculated.\n");
    if (thetafile_ != 0) mprintf("\tThetas output to %s\n", thetafile_->DataFilename().full());
  }

  init.DSL().SetDataSetsPending(true);
  masterDSL_ = init.DslPtr();
  return Action::OK;
}

// Action_MultiPucker::Setup()
Action::RetType Action_MultiPucker::Setup(ActionSetup& setup)
{
  Range actualRange;
  // If range is empty (i.e. no resrange arg given) look through all 
  // solute residues.
  if (resRange_.Empty())
    actualRange = setup.Top().SoluteResidues();
  else {
    // If user range specified, create new range shifted by -1 since internal
    // resnums start from 0.
    actualRange = resRange_;
    actualRange.ShiftBy(-1);
  }
  // Exit if no residues specified
  if (actualRange.Empty()) {
    mprinterr("Error: No residues specified for %s\n",setup.Top().c_str());
    return Action::ERR;
  }
  // Search for specified puckers in each residue in the range
  if (puckerSearch_.FindPuckers(setup.Top(), actualRange))
    return Action::SKIP;
  mprintf("\tResults of search in residue range [%s] for types", resRange_.RangeArg());
  puckerSearch_.PrintTypes();
  mprintf(", %u puckers found.\n", puckerSearch_.Npuckers());

  // Print selected puckers, set up DataSets
  data_.clear();
  amp_.clear();
  theta_.clear();
  puckerMethods_.clear();

  if (dsetname_.empty())
    dsetname_ = masterDSL_->GenerateDefaultName("MPUCKER");
  for (Pucker::PuckerSearch::mask_it pucker = puckerSearch_.begin();
                                     pucker != puckerSearch_.end(); ++pucker)
  {
    // setup/check the method
    Pucker::Method methodToUse = defaultMethod_;
    if (methodToUse == Pucker::UNSPECIFIED) {
      if (pucker->Natoms() > 5)
        methodToUse = Pucker::CREMER_POPLE;
      else
        methodToUse = Pucker::ALTONA_SUNDARALINGAM;
    } else if (methodToUse == Pucker::ALTONA_SUNDARALINGAM) {
      if (pucker->Natoms() > 5) {
        mprinterr("Error: Pucker '%s' res %i has too many atoms for Altona-Sundaralingam method.\n",
                  pucker->Name().c_str(), pucker->ResNum());
        return Action::ERR;
      }
    }
    puckerMethods_.push_back( methodToUse );

    int resNum = pucker->ResNum() + 1;
    // See if Dataset already present.
    MetaData md( dsetname_, pucker->Name(), resNum );
    DataSet* ds = masterDSL_->CheckForSet(md);
    if (ds == 0) {
      // Create new DataSet
      md.SetScalarMode( MetaData::M_PUCKER );
      md.SetScalarType( MetaData::PUCKER   ); // TODO pucker types
      ds = masterDSL_->AddSet( DataSet::DOUBLE, md );
      if (ds == 0) return Action::ERR;
      // Add to outfile
      if (outfile_ != 0)
        outfile_->AddDataSet( ds );
    }
    data_.push_back( ds );

    // Set up amplitude
    if (calc_amp_) {
      MetaData amp_md(dsetname_, pucker->Name() + "Amp", resNum);
      ds = masterDSL_->CheckForSet(amp_md);
      if (ds == 0) {
        md.SetScalarMode(MetaData::M_PUCKER);
        ds = masterDSL_->AddSet(DataSet::DOUBLE, amp_md);
        if (ds == 0) return Action::ERR;
        if (ampfile_ != 0) ampfile_->AddDataSet( ds );
      }
      amp_.push_back( ds );
    } else
      amp_.push_back( 0 );

    // Set up theta (> 5 atoms only)
    if (calc_theta_) {
      if (pucker->Natoms() < 6) {
        mprintf("Warning: 'theta' calc. not supported for < 6 atoms.\n");
        theta_.push_back( 0 );
      } else {
        MetaData theta_md(dsetname_, pucker->Name() + "Theta", resNum);
        ds = masterDSL_->CheckForSet(theta_md);
        if (ds == 0) {
          md.SetScalarMode(MetaData::M_PUCKER);
          ds = masterDSL_->AddSet(DataSet::DOUBLE, theta_md);
          if (ds == 0) return Action::ERR;
          if (thetafile_ != 0) thetafile_->AddDataSet( ds );
        }
        theta_.push_back( ds );
      }
    } else
      theta_.push_back( 0 );

    //if (debug_ > 0) {
      static const char* methodStr[] = { "Altona", "Cremer", "Unspecified" };
      mprintf("\t%zu [%s]: %s (%s)\n", data_.size(),ds->legend(), pucker->PuckerMaskString(setup.Top()).c_str(),
              methodStr[puckerMethods_.back()]);
    //}
  }
  return Action::OK;
}

// Action_MultiPucker::DoAction()
Action::RetType Action_MultiPucker::DoAction(int frameNum, ActionFrame& frm)
{
  const double* XYZ[6];
  for (int i = 0; i != 6; i++)
    XYZ[i] = 0;
  double pval, aval, tval;

  //std::vector<DataSet*>::const_iterator ds = data_.begin();
  //std::vector<Pucker::Method>::const_iterator method = puckerMethods_.begin();
  for (unsigned int idx = 0; idx != puckerSearch_.Npuckers(); idx++)
  {
    Pucker::PuckerMask const& pucker = puckerSearch_.FoundPucker(idx);
    // Since puckers are always at least 5 atoms, just reinit the 6th coord
    XYZ[5] = 0;
    // Get pucker coordinates
    unsigned int jdx = 0;
    for (Pucker::PuckerMask::atom_it atm = pucker.begin(); atm != pucker.end(); ++atm, ++jdx)
      XYZ[jdx] = frm.Frm().XYZ( *atm );
    // Do pucker calculation
    switch (puckerMethods_[idx]) {
      case Pucker::ALTONA_SUNDARALINGAM:
        pval = Pucker_AS( XYZ[0], XYZ[1], XYZ[2], XYZ[3], XYZ[4], aval );
        break;
      case Pucker::CREMER_POPLE:
        pval = Pucker_CP( XYZ[0], XYZ[1], XYZ[2], XYZ[3], XYZ[4], XYZ[5],
                          pucker.Natoms(), aval, tval );
        break;
      case Pucker::UNSPECIFIED : // Sanity check
        return Action::ERR;
    }

    if (amp_[idx] != 0) {
      aval *= Constants::RADDEG;
      amp_[idx]->Add(frameNum, &aval);
    }

    if (theta_[idx] != 0) {
      tval *= Constants::RADDEG;
      theta_[idx]->Add(frameNum, &tval);
    }

    pval = (pval * Constants::RADDEG) + offset_;
    if (pval > puckerMax_)
      pval -= PERIOD_;
    else if (pval < puckerMin_)
      pval += PERIOD_;
    data_[idx]->Add(frameNum, &pval);
  }
  return Action::OK;
}
