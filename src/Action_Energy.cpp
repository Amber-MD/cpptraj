#include "Action_Energy.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Energy::Action_Energy() :
  currentParm_(0)
{}

void Action_Energy::Help() {
  mprintf("\t[<name>] [<mask1>] [out <filename>]\n"
          "  Calculate energy for atoms in mask.\n");
}

static const char* Estring[] = {"bond", "angle", "dih", "vdw14", "elec14", "vdw", "elec", "total"};

// Action_Energy::Init()
Action::RetType Action_Energy::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  ENE_.SetDebug( debugIn );
  // Get keywords
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );

  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );

  // DataSet
  std::string setname = actionArgs.GetStringNext();
  if (setname.empty())
    setname = DSL->GenerateDefaultName("ENE");
  Energy_.clear(); 
  for (int i = 0; i <= (int)TOTAL; i++) {
    Energy_.push_back( DSL->AddSetAspect(DataSet::DOUBLE, setname, Estring[i]) );
    if (Energy_.back() == 0) return Action::ERR;
    // Add DataSet to DataFileList
    if (outfile != 0) outfile->AddSet( Energy_.back() );
  }

  mprintf("    ENERGY: Calculating energy for atoms in mask '%s'\n", Mask1_.MaskString());

  return Action::OK;
}

// Action_Energy::Setup()
/** Set angle up for this parmtop. Get masks etc.
  */
Action::RetType Action_Energy::Setup(Topology* currentParm, Topology** parmAddress) {
  if (currentParm->SetupCharMask(Mask1_)) return Action::ERR;
  if (Mask1_.None()) {
    mprinterr("Warning: Mask '%s' selects no atoms.\n", Mask1_.MaskString());
    return Action::ERR;
  }
  Mask1_.MaskInfo();
  Imask_ = Mask1_;
  Imask_.ConvertToIntMask();
  currentParm_ = currentParm;
  return Action::OK;
}

// Action_Energy::DoAction()
Action::RetType Action_Energy::DoAction(int frameNum, Frame* currentFrame,
                                            Frame** frameAddress)
{
  double Ebond = ENE_.E_bond(*currentFrame, *currentParm_, Mask1_);
  double Eangle = ENE_.E_angle(*currentFrame, *currentParm_, Mask1_);
  double Edih = ENE_.E_torsion(*currentFrame, *currentParm_, Mask1_);
  double Eq14 = 0.0;
  double Ev14 = ENE_.E_14_Nonbond(*currentFrame, *currentParm_, Mask1_, Eq14);
  double Eelec = 0.0;
  double Evdw = ENE_.E_Nonbond(*currentFrame, *currentParm_, Imask_, Eelec);
  double Etot = Ebond + Eangle + Edih + Eq14 + Ev14 + Evdw + Eelec;

  Energy_[BOND]->Add(frameNum, &Ebond);
  Energy_[ANGLE]->Add(frameNum, &Eangle);
  Energy_[DIHEDRAL]->Add(frameNum, &Edih);
  Energy_[V14]->Add(frameNum, &Ev14);
  Energy_[Q14]->Add(frameNum, &Eq14);
  Energy_[VDW]->Add(frameNum, &Evdw);
  Energy_[ELEC]->Add(frameNum, &Eelec);
  Energy_[TOTAL]->Add(frameNum, &Etot);

  return Action::OK;
}
