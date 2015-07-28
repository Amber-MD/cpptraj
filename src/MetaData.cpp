#include "MetaData.h"
#include "CpptrajStdio.h" // FIXME remove dependency
#include "StringRoutines.h"

const char* MetaData::Smodes[] = {"distance","angle","torsion","pucker","rms","matrix",0};
const char* MetaData::Stypes[] = {
  // Torsions
  "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "pucker",
  "chi",   "h1p",  "c2p",   "phi",   "psi",     "pchi", "omega",
  // Distance
  "noe",
  // Matrix
  "distance",    "covariance",          "mass-weighted covariance",
  "correlation", "distance covariance", "IDEA",
  "IRED",        "dihedral covariance",
  0 };
const MetaData::scalarMode MetaData::TypeModes[] = {
  M_TORSION, M_TORSION, M_TORSION, M_TORSION, M_TORSION, M_TORSION, M_PUCKER,
  M_TORSION, M_TORSION, M_TORSION, M_TORSION, M_TORSION, M_TORSION, M_TORSION,
  M_DISTANCE,
  M_MATRIX,  M_MATRIX, M_MATRIX,
  M_MATRIX,  M_MATRIX, M_MATRIX,
  M_MATRIX,  M_MATRIX,
  UNKNOWN_MODE };

std::string MetaData::ScalarDescription() const {
  std::string out("");
  if (scalarmode_ != UNKNOWN_MODE)
    out.append(", " + std::string(Smodes[scalarmode_]));
  if (scalartype_ != UNDEFINED)
    out.append("(" + std::string(Stypes[scalartype_]) + ")");
  return out;
}

MetaData::scalarMode MetaData::ModeFromKeyword(std::string const& key) {
  for (int i = 0; i != (int)UNKNOWN_MODE; i++)
    if (key.compare( Smodes[i] ) == 0) return (scalarMode)i;
  return UNKNOWN_MODE;
}

MetaData::scalarType MetaData::TypeFromKeyword(std::string const& key, scalarMode const& mIn)
{
  scalarMode dm = mIn;
  return TypeFromKeyword(key, dm);
}

MetaData::scalarType MetaData::TypeFromKeyword(std::string const& key, scalarMode& modeIn) {
  for (int i = 0; i != (int)UNDEFINED; i++)
    if (key.compare( Stypes[i] ) == 0) {
      if (modeIn != UNKNOWN_MODE) {
        // Is type valid for given mode?
        if (modeIn != TypeModes[i]) {
          mprinterr("Error: Type '%s' not valid for mode '%s'\n",Stypes[i],Smodes[TypeModes[i]]);
          return UNDEFINED;
        }
      } else
        modeIn = TypeModes[i];
      return (scalarType)i;
    }
  return UNDEFINED;
}

void MetaData::SetDefaultLegend() {
  if (!aspect_.empty() && idx_ == -1)
    legend_ = name_ + "[" + aspect_ + "]";
  else if (aspect_.empty() && idx_ != -1)
    legend_ = name_ + ":" + integerToString( idx_ );
  else if (!aspect_.empty() && idx_ != -1)
    legend_ = aspect_ + ":" + integerToString( idx_ );
  else
    legend_ = name_;
  if (ensembleNum_ != -1)
    legend_ += ("%" + integerToString( ensembleNum_ ));
}

std::string MetaData::PrintName() const {
  std::string out( name_ );
  if (!aspect_.empty())
    out.append("[" + aspect_ + "]");
  if (idx_ != -1)
    out.append(":" + integerToString(idx_));
  if (ensembleNum_ != -1)
    out.append("%" + integerToString(ensembleNum_));
  return out;
}

bool MetaData::Match_Exact(MetaData const& In) const {
  if (In.name_        != name_       ) return false;
  if (In.aspect_      != aspect_     ) return false;
  if (In.idx_         != idx_        ) return false;
  if (In.ensembleNum_ != ensembleNum_) return false;
  return true;
}

