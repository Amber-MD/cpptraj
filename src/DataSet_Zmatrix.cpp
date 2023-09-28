#include "DataSet_Zmatrix.h"
#include "CpptrajStdio.h"
#include "Structure/Zmatrix.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
DataSet_Zmatrix::DataSet_Zmatrix() : zmatrix_(0)
{

}

/** DESTRUCTOR */
DataSet_Zmatrix::~DataSet_Zmatrix() {
  if (zmatrix_ != 0) delete zmatrix_;
}
