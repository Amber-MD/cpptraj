#include "Analysis_Project.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"
#include "DataSet_Modes.h"
#include "StringRoutines.h"

// Analysis_Project::Help()
void Analysis_Project::Help() const {
  mprintf("\tevecs <evecs dataset> [name <name>] [out <outfile>] [beg <beg>] [end <end>]\n"
          "\t{[dihedrals <dataset arg>] | [data <dataset arg> ...]}\n"
          "  Calculate projection along given eigenvectors.\n");
}

// Analysis_Project::Setup()
Analysis::RetType Analysis_Project::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  beg_ = analyzeArgs.getKeyInt("beg", 1) - 1;
  end_ = analyzeArgs.getKeyInt("end", 2);

  std::string modesname = analyzeArgs.GetStringKey("evecs");
  if (modesname.empty()) {
    mprinterr("Error: No eigenvectors data set specified ('evecs <name>').\n");
    return Analysis::ERR;
  }

  // Check if DataSet exists
  modinfo_ = (DataSet_Modes*)setup.DSL().FindSetOfType( modesname, DataSet::MODES );
  if (modinfo_ == 0) {
    mprinterr("Error: No modes set '%s' found.\n", modesname.c_str());
    return Analysis::ERR;
  }

  // Output Filename
  DataFile* DF = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );

  Sets_.clear();
  std::string diharg = analyzeArgs.GetStringKey("dihedrals");
  std::string dataarg = analyzeArgs.GetStringKey("data");
  if (!diharg.empty() && !dataarg.empty()) {
    mprinterr("Error: Specify either 'dihedrals' or 'data', not both.\n");
    return Analysis::ERR;
  }
  if (!diharg.empty()) {
    // Get dihedral data sets
    Sets_.AddTorsionSets( setup.DSL().GetMultipleSets( diharg ) );
    if ( Sets_.empty() ) {
      mprinterr("Error: No valid data sets found.\n");
      return Analysis::ERR;
    }
  }
  if (!dataarg.empty()) {
    while (!dataarg.empty()) {
      if (Sets_.AppendSetsFromArgs( ArgList(dataarg), setup.DSL() )) {
        mprinterr("Error: Could not add data sets using argument '%s'\n", dataarg.c_str());
        return Analysis::ERR;
      }
      dataarg = analyzeArgs.GetStringKey("data");
    }
  }
  if (Sets_.empty()) {
    mprinterr("Error: No input data sets.\n");
    return Analysis::ERR;
  }

  // Set up data sets
  std::string setname = analyzeArgs.GetStringKey("name");
  if (setname.empty())
    setname = setup.DSL().GenerateDefaultName("Proj");
  for (int mode = beg_; mode < end_; ++mode) {
    int imode = mode + 1;
    DataSet* dout = setup.DSL().AddSet( DataSet::FLOAT, MetaData(setname, imode) );
    if (dout == 0) {
      mprinterr("Error: Could not create output dataset for mode %i\n", imode);
      return Analysis::ERR;
    }
    dout->SetLegend("Mode"+integerToString(imode));
    project_.push_back( dout );
    if (DF != 0) DF->AddDataSet( dout );
  }

  mprintf("    PROJECTION: Calculating projection using eigenvectors %i to %i of %s\n",
          beg_+1, end_, modinfo_->legend());
  if (DF != 0)
    mprintf("\tResults are written to %s\n", DF->DataFilename().full());
  mprintf("\t%zu input data sets.\n", Sets_.size());

  return Analysis::OK;
}

// Analysis_Project::Analyze()
Analysis::RetType Analysis_Project::Analyze() {
  if (modinfo_->Nmodes() < 1) {
    mprinterr("Error: modes set '%s' is empty.\n", modinfo_->legend());
    return Analysis::ERR;
  }
  // Check if beg and end are in bounds.
  if (end_ > modinfo_->Nmodes()) {
    mprintf("Warning: 'end' %i is greater than # evecs (%i); setting end to %i\n",
            end_, modinfo_->Nmodes(), modinfo_->Nmodes());
    end_ = modinfo_->Nmodes();
  }
  if (beg_ < 0 || beg_ >= end_) {
    mprinterr("Error: 'beg' %i out of bounds.\n", beg_+1);
    return Analysis::ERR;
  }
  // Check Modes type
  if (modinfo_->Meta().ScalarType() == MetaData::DIHCOVAR) {
    if ((int)Sets_.size() * 2 != modinfo_->VectorSize()) {
      mprinterr("Error: Number of dihedral data sets %zu does not correspond to"
                " number of eigenvectors %i\n", Sets_.size()*2, modinfo_->VectorSize());
      return Analysis::ERR;
    } else if ((int)Sets_.size() * 2 != modinfo_->NavgCrd()) {
      mprinterr("Error: Number of dihedral data sets %zu does not correspond to"
                " number of average elements %i\n", Sets_.size()*2, modinfo_->NavgCrd());
      return Analysis::ERR;
    }
  } else if (modinfo_->Meta().ScalarType() == MetaData::DATACOVAR) {
    if ((int)Sets_.size() != modinfo_->VectorSize()) {
      mprinterr("Error: Number of data sets %zu does not correspond to"
                " number of eigenvectors %i\n", Sets_.size(), modinfo_->VectorSize());
      return Analysis::ERR;
    } else if ((int)Sets_.size() != modinfo_->NavgCrd()) {
      mprinterr("Error: Number of data sets %zu does not correspond to"
                " number of average elements %i\n", Sets_.size(), modinfo_->NavgCrd());
      return Analysis::ERR;
    }
  } else {
    mprinterr("Error: Can only do data or dihedral modes currently.\n");
    return Analysis::ERR; // FIXME this check should be earlier
  }

  // Check that sets have same size
  unsigned int maxFrames = 0;
  if (!Sets_.empty()) {
    maxFrames = Sets_.Array().front()->Size();
    for (Array1D::const_iterator it = Sets_.begin(); it != Sets_.end(); ++it)
    {
      if ((*it)->Size() != maxFrames) {
        mprinterr("Error: Set '%s' does not have same size (%zu) as first set (%u)\n",
                  (*it)->legend(), (*it)->Size(), maxFrames);
        return Analysis::ERR;
      }
    }
  } 
  // Loop over frames
  for (unsigned int idx = 0; idx < maxFrames; idx++) {
    // Always start at first eigenvector element of first mode.
    const double* Vec = modinfo_->Eigenvector(beg_);
    if (modinfo_->Meta().ScalarType() == MetaData::DIHCOVAR ) {
      for (int mode = beg_; mode < end_; ++mode) {
        DataSet_Modes::AvgIt Avg = modinfo_->AvgBegin();
        double proj = 0.0;
        for (Array1D::const_iterator dih = Sets_.begin();
                                     dih != Sets_.end(); ++dih)
        {
          double theta = (*dih)->Dval( idx ) * Constants::DEGRAD;
          proj += (cos(theta) - *(Avg++)) * Vec[0];
          proj += (sin(theta) - *(Avg++)) * Vec[1];
          Vec += 2;
        }
        // TODO: Convert to degrees?
        float fproj = (float)proj;
        project_[mode]->Add( idx, &fproj );
      }
    } else if (modinfo_->Meta().ScalarType() == MetaData::DATACOVAR ) {
      for (int mode = beg_; mode < end_; ++mode) {
        DataSet_Modes::AvgIt Avg = modinfo_->AvgBegin();
        double proj = 0.0;
        for (Array1D::const_iterator it = Sets_.begin();
                                     it != Sets_.end(); ++it)
        {
          double dval = (*it)->Dval( idx ) - *(Avg++);
          proj += (dval * *Vec);
          Vec++;
        }
        float fproj = (float)proj;
        project_[mode]->Add( idx, &fproj );
      }
    }
  } // END loop over frames

  return Analysis::OK;
}
