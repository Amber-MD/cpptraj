#include "MetaData.h"
#include "CpptrajStdio.h" // FIXME remove dependency
#include "StringRoutines.h"

const char* MetaData::Smodes[] = {"distance","angle","torsion","pucker","rms","matrix","vector",0};
const char* MetaData::Stypes[] = {
  // Torsions
  "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "nu1",   "nu2",
  "h1p",   "c2p",  "chin",  "phi",   "psi",     "chip", "omega",
  "chi2",  "chi3", "chi4",  "chi5",
  // Pucker
  "pucker",
  // Distance
  "noe",
  // Matrix
  "distance",    "covariance",          "mass-weighted covariance",
  "correlation", "distance covariance", "IDEA",
  "IRED",        "dihedral covariance",
  // Vector
  "IRED",
  0 };
const MetaData::scalarMode MetaData::TypeModes[] = {
  M_TORSION, M_TORSION, M_TORSION, M_TORSION, M_TORSION, M_TORSION, M_TORSION, M_TORSION,
  M_TORSION, M_TORSION, M_TORSION, M_TORSION, M_TORSION, M_TORSION, M_TORSION,
  M_TORSION, M_TORSION, M_TORSION, M_TORSION,
  M_PUCKER,
  M_DISTANCE,
  M_MATRIX,  M_MATRIX, M_MATRIX,
  M_MATRIX,  M_MATRIX, M_MATRIX,
  M_MATRIX,  M_MATRIX,
  M_VECTOR,
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
  if (name_.empty() && !fileName_.empty())
    out.assign( fileName_.Full() );
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
  if (In.fileName_.Full() != fileName_.Full()) return false;
  if (In.fileName_.Base() != fileName_.Base()) return false;
  if (In.aspect_      != aspect_     ) return false;
  if (In.idx_         != idx_        ) return false;
  if (In.ensembleNum_ != ensembleNum_) return false;
  return true;
}

/** This version allows wildcards and ranges. */
bool MetaData::Match_WildCard(SearchString const& search) const
{
  //mprintf("DEBUG: Input: %s[%s]:%s  This Set: %s[%s]:%i\n",
  //        dsname.c_str(), aspect.c_str(), idxRange.RangeArg(), 
  //        name_.c_str(), aspect_.c_str(), idx_);
  // Match type if specified
  if (fileName_.empty()) {
    // No filename. Match name if specified.
    if ( WildcardMatch(search.NameArg(), name_) == 0 ) return false;
  } else {
    // Has associated file name. Exit if name specified no match with
    // name/base name/full path.
    if ( WildcardMatch(search.NameArg(), name_) == 0 &&
         !fileName_.MatchFullOrBase( search.NameArg() ) ) return false;
  }
  // If aspect specified make sure it matches.
  if ( WildcardMatch(search.AspectArg(), aspect_) == 0 ) return false;
  // Currently match any index if not specified.
  if (search.IdxRange().Front() != -1 && !search.IdxRange().InRange(idx_))
    return false;
  // Match any ensemble if not specified
  if (search.MemberRange().Front() != -1 && !search.MemberRange().InRange(ensembleNum_))
    return false;
  // If no aspect specified but dataset has aspect do not match.
  //if (aspect.empty() && !aspect_.empty()) return false;
  //mprintf("\tMATCH\n");
  return true;
}

// -----------------------------------------------------------------------------
/** Separate argument into metadata. Format: <name>[<aspect>]:<idx>%<ensemble>.
  * <name> and <aspect> support character wildcards * and ?; <idx> and
  * <ensemble> support ranges (e.g. 1-5,7,11).
  */
int MetaData::SearchString::ParseArgString(std::string const& argIn)
{
  name_arg_.assign( argIn );
  aspect_arg_.clear();
  idx_range_.Clear();
  member_range_.Clear();
  std::string idx_arg, member_arg;
  //mprinterr("DBG: ParseArgString called with %s\n", nameIn.c_str());
  // Separate out ensemble arg if present
  size_t idx_pos = name_arg_.find( '%' );
  if ( idx_pos != std::string::npos ) {
    // Advance to after the '.'
    member_arg = name_arg_.substr( idx_pos + 1 );
    // Drop the ensemble arg
    name_arg_.resize( idx_pos );
  }
  // Separate out index arg if present
  idx_pos = name_arg_.find( ':' );
  if ( idx_pos != std::string::npos ) {
    // Advance to after the ':'
    idx_arg = name_arg_.substr( idx_pos + 1 );
    //mprinterr("DBG:\t\tIndex Arg [%s]\n", idx_arg.c_str());
    // Drop the index arg
    name_arg_.resize( idx_pos );
  }
  // Separate out aspect arg if present
  size_t attr_pos0 = name_arg_.find_first_of( '[' );
  size_t attr_pos1 = name_arg_.find_last_of( ']' );
  if ( attr_pos0 != std::string::npos && attr_pos1 != std::string::npos ) {
    if ( (attr_pos0 != std::string::npos && attr_pos1 == std::string::npos) ||
         (attr_pos0 == std::string::npos && attr_pos1 != std::string::npos) )
    {
      mprinterr("Error: Malformed attribute ([<attr>]) in dataset arg %s\n", argIn.c_str());
      return 1;
    }
    // If '[' is at very beginning, this is a tag. Otherwise attribute.
    if (attr_pos0 > 0) {
      // Advance to after '[', length is position of ']' minus '[' minus 1 
      aspect_arg_ = name_arg_.substr( attr_pos0 + 1, attr_pos1 - attr_pos0 - 1 );
      //mprinterr("DBG:\t\tAttr Arg [%s]\n", aspect_arg_.c_str());
      // Drop the attribute arg
      name_arg_.resize( attr_pos0 );
    }
  }
  //mprinterr("DBG:\t\tName Arg [%s]\n", name_arg_.c_str());
  // If index arg is empty make wildcard (-1)
  if (idx_arg.empty() || idx_arg == "*")
    idx_range_.SetRange( -1, 0 );
  else
    idx_range_.SetRange( idx_arg );
  // If member arg is empty make none (-1)
  if (member_arg.empty() || member_arg == "*")
    member_range_.SetRange( -1, 0 );
  else
    member_range_.SetRange( member_arg );
  // If aspect arg not set and name is wildcard, make attribute wildcard.
  if (aspect_arg_.empty() && name_arg_ == "*")
    aspect_arg_.assign("*");

  return 0;
}
