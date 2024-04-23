#include "Analysis_Project.h"
#include "CpptrajStdio.h"
#include "DataSet_Modes.h"

// Analysis_Project::Help()
void Analysis_Project::Help() const {

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
  }

  return Analysis::OK;
}
