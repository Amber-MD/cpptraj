#include "Action_Vector.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Vector::Action_Vector() :
  Vector_(0)
{}

int Action_Vector::init() {
  // Keywords
  filename_ = actionArgs.GetStringKey("out");
 
  Vector_ = new VectorType();
  if (Vector_->Init( actionArgs )) {
    delete Vector_;
    return 1;
  }

  if (Vector_->Mode() == VectorType::VECTOR_NOOP) {
    mprinterr("Error: No vector mode specified.\n");
    return 1;
  }

  // Additional CORRIRED setup; check if modes need to be read or have been
  // previously read in by another VectorType.
  VectorType *Vtmp;
  if (Vector_->Mode() == VectorType::VECTOR_CORRIRED) {
    DSL->VectorBegin();
    while ( (Vtmp = DSL->NextVector()) != 0 ) {
      if ( *Vector_ == *Vtmp ) {
        if (Vector_->AssignMaster( Vtmp )) {
          mprinterr("Error: Could not assign vector master for CORRIRED.\n");
          delete Vector_;
          return 1;
        }
      }
    }
    // Load modes from file. modesfile name should be set by VectorType::Init.
    if (Vector_->NoModeInfo())
      if (Vector_->ReadModesFromFile()) {
        delete Vector_;
        return 1;
      }
  }

  // Allocate vector
  if ( Vector_->Allocate(DSL->MaxFrames()) ) return 1;

  // Add vector to datasetlist
  // TODO: Check for name conflicts
  DSL->AddDataSet( (DataSet*)Vector_ );

  Vector_->Info();

  // Check if output is supported for the current vector mode
  if (Vector_->Mode() == VectorType::VECTOR_CORRPLANE ||
      Vector_->Mode() == VectorType::VECTOR_CORR ||
      Vector_->Mode() == VectorType::VECTOR_CORRIRED ||
      Vector_->Mode() == VectorType::VECTOR_IRED)
  {
    if (!filename_.empty()) {
      mprintf("\tWarning: Output of corr, ired, corrired or corrplane vectors is not yet supported!\n");
      filename_.clear();
    }
  }

  if (!filename_.empty()) {
    mprintf("\tOutput will be dumped to a file, %s\n", filename_.c_str());
    DFL->Add( filename_.c_str(), (DataSet*)Vector_ );
  }

  return 0;
}

void Action_Vector::print() {
  //Vector_->Print();
}

int Action_Vector::setup() {
  return Vector_->Setup( currentParm );
}

int Action_Vector::action() {
  int err;

  switch (Vector_->Mode()) {
    case VectorType::VECTOR_CORRIRED:
    case VectorType::VECTOR_CORR:
    case VectorType::VECTOR_CORRPLANE:
      err = Vector_->Action_CORR( currentFrame );
      break;
    case VectorType::VECTOR_DIPOLE:
      err = Vector_->Action_DIPOLE( currentFrame, currentParm );
      break;
    case VectorType::VECTOR_PRINCIPAL_X:
    case VectorType::VECTOR_PRINCIPAL_Y:
    case VectorType::VECTOR_PRINCIPAL_Z:
      err = Vector_->Action_PRINCIPAL( currentFrame ); 
      break;
    case VectorType::VECTOR_MASK:
      err = Vector_->Action_MASK( currentFrame );
      break;
    case VectorType::VECTOR_IRED:
      err = Vector_->Action_IRED( currentFrame );
      break;
    case VectorType::VECTOR_BOX:
      err = Vector_->Action_BOX( currentFrame );
      break;
    default:
      mprinterr("Error: Unhandled vector operation.\n");
      return 1;
  }
  return err;
}

