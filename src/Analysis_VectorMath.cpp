#include <algorithm> // std::min, std::max
#include "Analysis_VectorMath.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "DataSet_double.h"

/// Strings corresponding to modes, used in output.
const char* Analysis_VectorMath::ModeString[] = {
  "Dot product", "Angle from dot product", "Cross product" };

// CONSTRUCTOR
Analysis_VectorMath::Analysis_VectorMath() :
  mode_(DOTPRODUCT), vinfo1_(0), vinfo2_(0), DataOut_(0), norm_(false) {}

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
  vinfo1_ = (DataSet_Vector*)setup.DSL().FindSetOfType( analyzeArgs.GetStringKey("vec1"),
                                                   DataSet::VECTOR );
  vinfo2_ = (DataSet_Vector*)setup.DSL().FindSetOfType( analyzeArgs.GetStringKey("vec2"),
                                                   DataSet::VECTOR );
  if (vinfo1_ == 0 ) {
    mprinterr("Error: 'vec1' not found.\n");
    return Analysis::ERR;
  }
  if (vinfo2_ == 0) {
    mprinterr("Error: 'vec2' not found.\n");
    return Analysis::ERR;
  }
  std::string setname = analyzeArgs.GetStringKey("name");
  norm_ = analyzeArgs.hasKey("norm");
  // Check for dotproduct/crossproduct keywords. Default is dotproduct.
  DataOut_ = 0;
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
  // Set up output data set based on mode
  DataOut_ = setup.DSL().AddSet(dtype, setname, dname);
  if (DataOut_ == 0) return Analysis::ERR;
  // Set up output file in DataFileList if necessary
  DataFile* outfile = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  if (outfile != 0) outfile->AddDataSet( DataOut_ );

  // Print Status
  mprintf("    VECTORMATH: Calculating %s of vectors %s and %s\n", 
            ModeString[mode_], vinfo1_->legend(), vinfo2_->legend());
  if (norm_) mprintf("\tVectors will be normalized.\n");
  if (outfile != 0)
    mprintf("\tResults are written to %s\n", outfile->DataFilename().full());

  return Analysis::OK;
}

// Analysis_VectorMath::DotProduct()
int Analysis_VectorMath::DotProduct(unsigned int vmax, unsigned int v1inc, unsigned int v2inc)
{
  DataSet_double& Out = static_cast<DataSet_double&>( *DataOut_ );
  DataSet_Vector& V1 = static_cast<DataSet_Vector&>( *vinfo1_ );
  DataSet_Vector& V2 = static_cast<DataSet_Vector&>( *vinfo2_ );
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
int Analysis_VectorMath::CrossProduct(unsigned int vmax, unsigned int v1inc, unsigned int v2inc)
{
  DataSet_Vector& Out = static_cast<DataSet_Vector&>( *DataOut_ );
  DataSet_Vector& V1 = static_cast<DataSet_Vector&>( *vinfo1_ );
  DataSet_Vector& V2 = static_cast<DataSet_Vector&>( *vinfo2_ );
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

// Analysis_VectorMath::Analyze()
Analysis::RetType Analysis_VectorMath::Analyze() {
  if (vinfo1_->Size() == 0 || vinfo2_->Size() == 0) {
    mprinterr("Error: One or both vectors is empty.\n");
    return Analysis::ERR;
  }
  // Either vec1 or vec2 can be size 1, else both need to have same size.
  unsigned int v1inc = 1;
  unsigned int v2inc = 1;
  if (vinfo1_->Size() != vinfo2_->Size()) {
    if (vinfo1_->Size() == 1)
      v1inc = 0;
    else if (vinfo2_->Size() == 1)
      v2inc = 0;
    else {
      mprinterr("Error: # Frames in vec %s (%i) != # Frames in vec %s (%i)\n",
                vinfo1_->legend(), vinfo1_->Size(),
                vinfo2_->legend(), vinfo2_->Size());
      return Analysis::ERR;
    }
  }
  unsigned int vmax = std::max( vinfo1_->Size(), vinfo2_->Size() );
  mprintf("\t'%s' size %zu, '%s' size %zu, output size %zu\n",
          vinfo1_->legend(), vinfo1_->Size(), vinfo2_->legend(), vinfo2_->Size(), vmax);
  int err = 0;
  if (mode_ == CROSSPRODUCT)
    err = CrossProduct(vmax, v1inc, v2inc);
  else // DOTPRODUCT || DOTANGLE
    err = DotProduct(vmax, v1inc, v2inc);
  if (err != 0) return Analysis::ERR;
  return Analysis::OK;
}
