#include "Analysis_CrankShaft.h"
#include "CpptrajStdio.h"

Analysis_CrankShaft::Analysis_CrankShaft() :
  start_(0),
  stop_(-1),
  offset_(1),
  type_(ANGLE),
  scalar1_(0),
  scalar2_(0)
{}

/** analyze crankshaft {angle | distance} <scalar-name1> <scalar-name2> 
  *                     info <string>
  */
int Analysis_CrankShaft::Setup(DataSetList *DSLin) {
  info_ = analyzeArgs_.GetStringKey("info");
  if (info_.empty())
    info_.assign("");

  if (analyzeArgs_.hasKey("angle"))
    type_ = ANGLE;
  else if (analyzeArgs_.hasKey("distance"))
    type_ = DISTANCE;

  filename_ = analyzeArgs_.GetStringKey("out");

  start_ = analyzeArgs_.getKeyInt("start", 1);
  --start_;
  stop_ = analyzeArgs_.getKeyInt("stop", -1);
  offset_ = analyzeArgs_.getKeyInt("offset",1);

  // Get dataset names
  std::string name1 = analyzeArgs_.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: crankshaft: No name specified for dataset 1.\n");
    return 1;
  }
  std::string name2 = analyzeArgs_.GetStringNext();
  if (name2.empty()) {
    mprinterr("Error: crankshaft: No name specified for dataset 2.\n");
    return 1;
  }

  // Get datasets
  scalar1_ = DSLin->Get( name1.c_str() );
  if (scalar1_ == NULL) {
    mprinterr("Error: crankshaft: Dataset %s not found.\n", name1.c_str());
    return 1;
  }
  scalar2_ = DSLin->Get( name2.c_str() );
  if (scalar2_ == NULL) {
    mprinterr("Error: crankshaft: Dataset %s not found.\n", name2.c_str());
    return 1;
  } 

  // INFO:
  mprintf("    ANALYZE CRANKSHAFT: %s ", info_.c_str());
  if (type_==ANGLE)
    mprintf("angles");
  else if (type_==DISTANCE)
    mprintf("distances");
  mprintf(" named %s and %s\n", name1.c_str(), name2.c_str());
  mprintf("\tFrames %i to ", start_+1);
  if (stop_==-1)
    mprintf("last");
  else
    mprintf("%i", stop_);
  mprintf(", offset %i\n", offset_);

  return 0;
}

int Analysis_CrankShaft::Analyze() {
  double v1_avg[6][6], v2_avg[6][6], v1_sd[6][6], v2_sd[6][6];
  int visits[6][6], transitions[6][6];

  // Check that scalar1 and scalar2 have same # data points.
  int Nelements = scalar1_->Size();
  if (Nelements != scalar2_->Size()) {
    mprinterr("Error: crankshaft: # elements in dataset %s (%i) not equal to\n",
              scalar1_->c_str(), Nelements);
    mprinterr("                   # elements in dataset %s (%i)\n",
              scalar2_->c_str(), scalar2_->Size());
    return 1;
  }
  if (stop_ == -1)
    stop_ = Nelements;
  if (start_ >= Nelements) {
    mprinterr("Error: crankshaft: start (%i) >= total # elements.\n",start_+1, Nelements);
    return 1;
  }
  int totalFrames = (stop_ - start_) / offset_;
  mprintf("\tcrankshaft: Processing %i frames.\n", totalFrames);

  // Initialization
  for (int j=0;j<6;j++) {
    for (int k=0;k<6;k++) {
      v1_avg[j][k] = 0.0;
      v2_avg[j][k] = 0.0;
      v1_sd[j][k] = 0.0;
      v2_sd[j][k] = 0.0;

      visits[j][k] = 0;
      transitions[j][k] = 0;
    }
  }

  // Open output file
  CpptrajFile outfile;
  if (outfile.OpenWrite( filename_ )) return 1;

  // MAIN LOOP over frames
  int i1 = 0;
  int i2 = 0;
  int initial_i1 = 0;
  int initial_i2 = 0;
  int final_i1 = 0;
  int final_i2 = 0;
  double initial_v1 = 0;
  double initial_v2 = 0;
  double final_v1 = 0;
  double final_v2 = 0;
  for (int frame = start_; frame < stop_; frame += offset_) {
    double v1 = scalar1_->Dval( frame );
    double v2 = scalar1_->Dval( frame );

    if ( type_ == DISTANCE ) {
      // This is a DISTANCE
      i1 = (v1 - 1.0) / 1;     // The current algorithm is a test and aims to bin things 
      if (i1 > 5) i1 = 5;      // from 0->5 starting from values < 2A, in increments of 1
      i2 = (v2 - 1.0) / 1;     // angstrom to > 6 A.  i.e. value-1.0/1                   
      if (i2 > 5) i2 = 5;
    } else { // if type_ == ANGLE
      // This is an ANGLE
      //   -- bin from 0->5
      //   -- subtract 30 from value, such that 0->60 = g+, 60->120 = a+, 120->180 = t

      v1 -= 30.0;
      v2 -= 30.0;

      if (v1 < 0)    v1 += 360.0;
      if (v1 > 360)  v1 -= 360.0;
      if (v2 < 0)    v2 += 360.0;
      if (v2 > 360)  v2 -= 360.0;

      i1 = v1 / 60;
      i2 = v2 / 60;
    }

    // Store initial and final bins/values
    if (frame == start_) {
      initial_i1 = i1;
      initial_i2 = i2;
      initial_v1 = scalar1_->Dval( frame );
      initial_v2 = scalar2_->Dval( frame );
    }
    final_i1 = i1;
    final_i2 = i2;
    final_v1 = scalar1_->Dval( frame );
    final_v2 = scalar2_->Dval( frame );

  } // END loop over frames

  return 0;
}

