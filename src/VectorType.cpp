#include "VectorType.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
VectorType::VectorType() :
  totalFrames_(-1),
  frame_(0),
  mode_(VECTOR_NOOP),
  cx_(0),
  cy_(0),
  cz_(0),
  vx_(0),
  vy_(0),
  vz_(0),
  master_(0),
  ibeg_(1),
  iend_(50),
  order_(2),
  npair_(-1),
  avgcrd_(0),
  rave_(0),
  r3iave_(0),
  r6iave_(0),
  cftmp_(0),
  p2cftmp_(0),
  rcftmp_(0)
{}

int VectorType::Init(ArgList& argIn) {
  filename_ = argIn.GetStringKey("out");

  order_ = argIn.getKeyInt("order",2);
  if (order_ < 0 || order_ > 2) {
    mprintf("Warning: vector order out of bounds (<0 or >2), resetting to 2.\n");
    order_ = 2;
  }

  // Acceptable args: principal, principal x, principal y, principal z
  int pos = argIn.KeyPosition("principal");
  if (pos == argIn.Nargs() - 1) {
    argIn.MarkArg(pos);
    mode_ = VECTOR_PRINCIPAL_X;
  } else if (pos!=-1) {
    argIn.MarkArg(pos);
    mode_ = VECTOR_PRINCIPAL_X;
    char vecchar = argIn[pos+1][0];
    if (vecchar == 'x' ||
        vecchar == 'y' ||
        vecchar == 'z')
    {
      argIn.MarkArg(pos+1);
      switch (vecchar) {
        case 'x': mode_ = VECTOR_PRINCIPAL_X; break;
        case 'y': mode_ = VECTOR_PRINCIPAL_Y; break;
        case 'z': mode_ = VECTOR_PRINCIPAL_Z; break;
      }
    }
  } else if (argIn.hasKey("dipole"))
    mode_ = VECTOR_DIPOLE;
  else if (argIn.hasKey("box"))
    mode_ = VECTOR_BOX;
  else if (argIn.hasKey("corrplane"))
    mode_ = VECTOR_CORRPLANE;
  else if (argIn.hasKey("corrired"))
    mode_ = VECTOR_CORRIRED;
  else if (argIn.hasKey("corr"))
    mode_ = VECTOR_CORR;
  else if (argIn.hasKey("ired"))
    mode_ = VECTOR_IRED;
  else
    mode_ = VECTOR_MASK;

  // VECTOR_CORRIRED
  if (mode_ == VECTOR_CORRIRED) {
    // Get Pair number
    npair_ = argIn.getKeyInt("npair",0);
    if (npair_ == 0) {
      mprinterr("Error: VectorType::Init(): No 'npair <#>' arg given, needed for 'corrired'.\n");
      return 1;
    }
    char* modesfile = argIn.getKeyString("modes",NULL);
    if (modesfile==NULL) {
      mprinterr("Error: VectorType::Init(): No 'modes' <file> arg given, needed for 'corrired'.\n");
      return 1;
    }
    ibeg_ = argIn.getKeyInt("beg",1);
    iend_ = argIn.getKeyInt("end", 50);
    
  }

  return 0;
}
    
  
