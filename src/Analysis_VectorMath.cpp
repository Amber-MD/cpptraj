#include <algorithm> // std::min, std::max
#include "Analysis_VectorMath.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "DataSet_Vector.h"
#include "DataSet_double.h"

/// Strings corresponding to modes, used in output.
const char* Analysis_VectorMath::ModeString[] = {
  "Dot product", "Angle from dot product", "Cross product" };

// CONSTRUCTOR
Analysis_VectorMath::Analysis_VectorMath() :
  mode_(DOTPRODUCT),
  norm_(false)
{}

void Analysis_VectorMath::Help() const {
  mprintf("\tvec1 <vecname1> vec2 <vecname2> [out <filename>] [norm] [name <setname>]\n"
          "\t[ dotproduct | dotangle | crossproduct ]\n"
          "  Calculate dot product, angle from dot product (degrees), or cross product\n"
          "  for specified vectors. Either vec1 or vec2 can be size 1, otherwise they\n"
          "  must both be the same size.\n");
}

// Analysis_VectorMath::Setup()
Analysis::RetType Analysis_VectorMath::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  // Get Vectors
  DataSetList vsets1 = setup.DSL().GetSetsOfType( analyzeArgs.GetStringKey("vec1"),
                                                  DataSet::VECTOR );
  if (vsets1.empty()) {
    mprinterr("Error: 'vec1' not found.\n");
    return Analysis::ERR;
  }
  DataSetList vsets2 = setup.DSL().GetSetsOfType( analyzeArgs.GetStringKey("vec2"),
                                                  DataSet::VECTOR );
  if (vsets2.empty()) {
    mprinterr("Error: 'vec2' not found.\n");
    return Analysis::ERR;
  }

  if (vsets1.size() != vsets2.size()) {
    mprinterr("Error: 'vec1' (%zu) and 'vec2' (%zu) do not select the same number of sets.\n",
              vsets1.size(), vsets2.size());
    return Analysis::ERR;
  }

  for (DataSetList::const_iterator it = vsets1.begin(); it != vsets1.end(); ++it)
    vinfo1_.push_back( static_cast<DataSet_Vector*>( *it ) );
  for (DataSetList::const_iterator it = vsets2.begin(); it != vsets2.end(); ++it)
    vinfo2_.push_back( static_cast<DataSet_Vector*>( *it ) );

  std::string setname = analyzeArgs.GetStringKey("name");
  norm_ = analyzeArgs.hasKey("norm");
  // Check for dotproduct/crossproduct keywords. Default is dotproduct.
  mode_ = DOTPRODUCT;
  DataSet::DataType dtype = DataSet::DOUBLE;
  const char* dname = "Dot";
  if (analyzeArgs.hasKey("dotproduct")) {
    mode_ = DOTPRODUCT;
  } else if (analyzeArgs.hasKey("dotangle")) {
    mode_ = DOTANGLE;
    norm_ = true; // Vecs must be normalized for angle calc to work
    dname = "Angle";
  } else if (analyzeArgs.hasKey("crossproduct")) {
    mode_ = CROSSPRODUCT;
    dtype = DataSet::VECTOR;
    dname = "Cross";
  }
  // Set up output file in DataFileList if necessary
  DataFile* outfile = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  // Set up output data sets based on mode
  if (setname.empty())
    setname = setup.DSL().GenerateDefaultName( dname );
  MetaData md(setname);
  int idx = -1;
  if (vinfo1_.size() > 1)
    idx = 0;
  for (unsigned int ii = 0; ii < vinfo1_.size(); ii++)
  {
    if (idx > -1)
      md.SetIdx( idx++ );
    DataSet* dsout = setup.DSL().AddSet(dtype, md);
    if (dsout == 0) return Analysis::ERR;
    if (outfile != 0) outfile->AddDataSet( dsout );
    DataOut_.push_back( dsout );
  }

  // Print Status
  mprintf("    VECTORMATH:");
  if (vinfo1_.size() == 1)
    mprintf(" Calculating %s of vectors %s and %s\n", ModeString[mode_], vinfo1_[0]->legend(), vinfo2_[0]->legend());
  else {
    mprintf(" Calculating %s of:\n", ModeString[mode_]);
    for (unsigned int ii = 0; ii < vinfo1_.size(); ii++)
      mprintf("\t  %s and %s\n", vinfo1_[ii]->legend(), vinfo2_[ii]->legend());
  }
  if (norm_) mprintf("\tVectors will be normalized.\n");
  if (outfile != 0)
    mprintf("\tResults are written to %s\n", outfile->DataFilename().full());

  return Analysis::OK;
}

// Analysis_VectorMath::DotProduct()
int Analysis_VectorMath::DotProduct(DataSet* Dout, DataSet_Vector& V1, DataSet_Vector& V2,
                                    unsigned int vmax, unsigned int v1inc, unsigned int v2inc)
const
{
  DataSet_double& Out = static_cast<DataSet_double&>( *Dout );
  Out.Resize( vmax );
  unsigned int v1 = 0;
  unsigned int v2 = 0;
  for (unsigned int v = 0; v < vmax; ++v, v1 += v1inc, v2 += v2inc) {
    if (norm_) {
      V1[v1].Normalize();
      V2[v2].Normalize();
    }
    if (mode_ == DOTPRODUCT)
      Out[v] = V1[v1] * V2[v2];
    else // DOTANGLE
      Out[v] = V1[v1].Angle( V2[v2] ) * Constants::RADDEG;
  }
  return 0;
}

// Analysis_VectorMath::CrossProduct()
int Analysis_VectorMath::CrossProduct(DataSet* Dout, DataSet_Vector& V1, DataSet_Vector& V2,
                                      unsigned int vmax, unsigned int v1inc, unsigned int v2inc)
const
{
  DataSet_Vector& Out = static_cast<DataSet_Vector&>( *Dout );
  Out.ReserveVecs( V1.Size() );
  unsigned int v1 = 0;
  unsigned int v2 = 0;
  for (unsigned int v = 0; v < vmax; ++v, v1 += v1inc, v2 += v2inc) {
    if (norm_) {
      V1[v1].Normalize();
      V2[v2].Normalize();
    }
    Out.AddVxyz( V1[v1].Cross( V2[v2] ) );
  }
  return 0;
}

/** Perform specified operation on two vectors, store result.
  * NOTE: Not passed in as const ref in case normalization is specified.
  */
int Analysis_VectorMath::DoMath(DataSet* Dout, DataSet_Vector& v1, DataSet_Vector& v2)
const
{
  if (v1.Size() == 0 || v2.Size() == 0) {
    mprinterr("Error: One or both vectors is empty.\n");
    return 1;
  }
  // Either vec1 or vec2 can be size 1, else both need to have same size.
  unsigned int v1inc = 1;
  unsigned int v2inc = 1;
  if (v1.Size() != v2.Size()) {
    if (v1.Size() == 1)
      v1inc = 0;
    else if (v2.Size() == 1)
      v2inc = 0;
    else {
      mprinterr("Error: # Frames in vec %s (%zu) != # Frames in vec %s (%zu)\n",
                v1.legend(), v1.Size(),
                v2.legend(), v2.Size());
      return 1;
    }
  }
  unsigned int vmax = std::max( v1.Size(), v2.Size() );
  mprintf("\t'%s' size %zu, '%s' size %zu, output size %u\n",
          v1.legend(), v1.Size(), v2.legend(), v2.Size(), vmax);
  int err = 0;
  if (mode_ == CROSSPRODUCT)
    err = CrossProduct(Dout, v1, v2, vmax, v1inc, v2inc);
  else // DOTPRODUCT || DOTANGLE
    err = DotProduct(Dout, v1, v2, vmax, v1inc, v2inc);
  return err;
}

// Analysis_VectorMath::Analyze()
Analysis::RetType Analysis_VectorMath::Analyze() {
  for (unsigned int ii = 0; ii < vinfo1_.size(); ii++)
  {
    int err = DoMath( DataOut_[ii], *(vinfo1_[ii]), *(vinfo2_[ii]) );
    if (err != 0) return Analysis::ERR;
  }
  return Analysis::OK;
}
