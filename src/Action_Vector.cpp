#include "Action_Vector.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Vector::Action_Vector() :
  Vector_(0)
{}

int Action_Vector::init() {
  
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

  return 0;
}

void Action_Vector::print() {
  Vector_->Print();
}

int Action_Vector::setup() {
  return Vector_->Setup( currentParm );
}

