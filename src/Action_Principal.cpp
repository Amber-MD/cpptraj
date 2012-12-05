#include "Action_Principal.h"
#include "CpptrajStdio.h"
#include "Matrix_3x3.h"

// CONSTRUCTOR
Action_Principal::Action_Principal() :
  doRotation_(false),
  useMass_(false),
  debug_(0)
{ }

void Action_Principal::Help() {
  mprintf("principal [<mask>] [dorotation] [mass]\n");
}

// Action_Principal::init()
Action::RetType Action_Principal::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Keywords
  doRotation_ = actionArgs.hasKey("dorotation");
  useMass_ = actionArgs.hasKey("mass");

  // Masks
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    PRINCIPAL");
  if (doRotation_)
    mprintf(" with rotation by");
  else
    mprintf(" without rotation by");
  if (useMass_)
    mprintf(" center of mass");
  else
    mprintf(" center of geometry");
  mprintf(", atoms selected by [%s]\n", mask_.MaskString());

  return Action::OK;
}

// Action_Principal::setup()
Action::RetType Action_Principal::Setup(Topology* currentParm, Topology** parmAddress) {
  if (currentParm->SetupIntegerMask(mask_)) return Action::ERR;

  if (mask_.None()) {
    mprintf("Warning: No atoms selected for %s [%s].\n",currentParm->c_str(), mask_.MaskString());
    return Action::ERR;
  }

  mprintf("\tSelected %i atoms.\n", mask_.Nselected());
  return Action::OK;
}

// Action_Principal::action()
Action::RetType Action_Principal::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Matrix_3x3 Inertia;
  Vec3 Eval;

  currentFrame->CalculateInertia( mask_, Inertia );
  //printMatrix_3x3("PRINCIPAL INERTIA", Inertia);

  // NOTE: Diagonalize_Sort_Chirality places sorted eigenvectors in rows.
  Inertia.Diagonalize_Sort_Chirality( Eval, debug_ );
  if (debug_ > 2) {
    Eval.Print("PRINCIPAL EIGENVALUES");
    //TEMP.Print("GENERAL");
    Inertia.Print("PRINCIPAL EIGENVECTORS (Rows)");
  }
  
  // Rotate - since Evec is already transposed (eigenvectors
  // are returned in rows) just do plain rotation to affect an
  // inverse rotation.
  if (doRotation_)
    currentFrame->Rotate( Inertia );

  return Action::OK;
}
