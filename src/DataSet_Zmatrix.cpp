#include "DataSet_Zmatrix.h"
#include "CpptrajStdio.h"
#include "Structure/Zmatrix.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
DataSet_Zmatrix::DataSet_Zmatrix() :
  DataSet(ZMATRIX, GENERIC, TextFormat(TextFormat::DOUBLE, 8, 3, 7), 1), // i j k l dist theta phi
  zmatrix_(0)
{
  zmatrix_ = new Zmatrix();
}

/** DESTRUCTOR */
DataSet_Zmatrix::~DataSet_Zmatrix() {
  if (zmatrix_ != 0) delete zmatrix_;
}

/** \return memory usage in bytes. */
size_t DataSet_Zmatrix::MemUsageInBytes() const {
  return zmatrix_->sizeInBytes();
}

/** \return number of internal coords */
size_t DataSet_Zmatrix::Size() const {
  return zmatrix_->N_IC();
}

/** Reserve space in the IC_ array. */
int DataSet_Zmatrix::Allocate( SizeArray const& sizeIn ) {
  if (!sizeIn.empty())
    zmatrix_->reserve( sizeIn[0] );
  return 0;
}

/** Write to file. */
void DataSet_Zmatrix::WriteBuffer(CpptrajFile &cbuffer, SizeArray const& frame) const {
  if (frame[0] >= zmatrix_->N_IC())
    cbuffer.Printf(format_.fmt(), 0, 0, 0, 0, 0, 0, 0);
  else {
    //mprintf("DEBUG: Format: %s\n", format_.fmt());
    InternalCoords const& ic = (*zmatrix_)[frame[0]];
    cbuffer.Printf(format_.fmt(),
                   (double)(ic.AtI()+1), (double)(ic.AtJ()+1),
                   (double)(ic.AtK()+1), (double)(ic.AtL()+1),
                   ic.Dist(), ic.Theta(), ic.Phi());
  }
}
