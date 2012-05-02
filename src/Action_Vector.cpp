#include "Action_Vector.h"

// CONSTRUCTOR
Action_Vector::Action_Vector() :
  Vector_(0)
{}

int Action_Vector::init() {
  
  Vector_ = new VectorType();
  Vector_->Init( actionArgs );

  // Additional CORRIRED setup; check if modes need to be read or have been
  // previously read in by another VectorType.
  if (Vector_->Type() == VectorType::VECTOR_CORRIRED) {
    DSL->VectorBegin();
    while ( (VectorType *Vtmp = DSL->NextVector()) != NULL ) {
      if ( *Vector_ == *Vtmp )
        Vector_->AssignModes( Vtmp );

  return 0;
}
