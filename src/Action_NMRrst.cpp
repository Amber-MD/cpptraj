#include <cmath>
#include "Action_NMRrst.h"
#include "DataSet_double.h"
#include "StringRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_NMRrst::Action_NMRrst() : useMass_(false), resOffset_(0) {} 

void Action_NMRrst::Help() {
  mprintf("\t[<name>] file <rstfile> [name <dataname>] [geom] [noimage] [resoffset <r>]\n");
}

/// \return true if character is not a 'skippable' one.
static inline bool NotSkipChar( const char* ptr ) {
  return (ptr != 0 && *ptr != '#' && *ptr != '!' && *ptr != '\n' && *ptr != '\r');
}

// Action_NMRrst::Init()
Action::RetType Action_NMRrst::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get Keywords
  Image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  useMass_ = !(actionArgs.hasKey("geom"));
  resOffset_ = actionArgs.getKeyInt("resoffset", 0);
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  std::string rstfilename = actionArgs.GetStringKey("file");
  if (rstfilename.empty()) {
    mprinterr("Error: must specify an NMR restraint filename with 'file <rstfile>'\n");
    return Action::ERR;
  }
  std::string setname = actionArgs.GetStringKey("name");
  if (setname.empty())
    setname = DSL->GenerateDefaultName("NMR");

  // Read in NMR restraints.
  BufferedLine infile;
  if (infile.OpenFileRead( rstfilename )) return Action::ERR;
  // Try to determine what kind of file.
  const char* ptr = infile.Line();
  // Try to skip past any blank lines and comments
  while ( NotSkipChar( ptr ) )
    ptr = infile.Line();
  if (ptr == 0) {
    mprinterr("Error: Unexpected end of restraint file.\n");
    return Action::ERR;
  }
  std::string inputLine( ptr );
  infile.CloseFile();
  // Re-open file
  if (infile.OpenFileRead( rstfilename )) return Action::ERR;
  int err = 0;
  if ( inputLine.compare(0, 7, "*HEADER") ||
       inputLine.compare(0, 6, "*TITLE") ||
       inputLine.compare(0, 6, "assign") )
    // XPLOR
    err = ReadXplor( infile );
  else
    // Assume DIANA/Amber
    err = ReadAmber( infile );
  if (err != 0) {
    mprinterr("Error: Could not parse restraint file.\n");
    return Action::ERR;
  }

  // Set up distances.
  int num_noe = 1;
  for (noeArray::iterator noe = NOEs_.begin(); noe != NOEs_.end(); ++noe, ++num_noe) {
     // Dataset to store distances
    (*noe).dist_ = DSL->AddSetIdxAspect(DataSet::DOUBLE, setname, num_noe, "NOE");
    if ((*noe).dist_==0) return Action::ERR;
    (*noe).dist_->SetScalar( DataSet::M_DISTANCE, DataSet::NOE );
    ((DataSet_double*)(*noe).dist_)->SetNOE((*noe).bound_, (*noe).boundh_, (*noe).rexp_);
    (*noe).dist_->SetLegend((*noe).dMask1_.MaskExpression() + " and " + (*noe).dMask2_.MaskExpression());
    // Add dataset to data file
    if (outfile != 0) outfile->AddSet( (*noe).dist_ );
  }
 
  mprintf("    NMRRST: %zu NOEs\n", NOEs_.size());
  mprintf("\tShifting residue numbers in restraint file by %i\n", resOffset_);
  // DEBUG - print NOEs
  for (noeArray::const_iterator noe = NOEs_.begin(); noe != NOEs_.end(); ++noe)
    mprintf("\t'%s' to '%s'  %f < %f < %f\n", (*noe).dMask1_.MaskString(),
            (*noe).dMask2_.MaskString(), (*noe).bound_, (*noe).rexp_, (*noe).boundh_);
  if (!Image_.UseImage()) 
    mprintf("\tNon-imaged");
  else
    mprintf("\tImaged");
  if (useMass_) 
    mprintf(", center of mass.\n");
  else
    mprintf(", geometric center.\n");

  return Action::OK;
}
// -----------------------------------------------------------------------------
int Action_NMRrst::ReadAmber( BufferedLine& infile ) {
  return 1;
}

/// Replace any '#' with '='
// FIXME: Is this legitimate?
inline void ReplaceChars( std::string& maskExpression ) {
  for (std::string::iterator c = maskExpression.begin();
                             c != maskExpression.end(); ++c)
    if (*c == '#') *c = '=';
}

static inline std::string GetSelection(ArgList& line, int offset) {
  std::string maskExpression = ":";
  int resnum = line.getKeyInt("resid",0) + offset;
  if (resnum < 1) return std::string("");
  maskExpression += integerToString(resnum);
  maskExpression += "@";
  maskExpression += line.GetStringKey("name");
  ReplaceChars( maskExpression );
  return maskExpression;
}

int Action_NMRrst::ReadXplor( BufferedLine& infile ) {
  noeDataType NOE;
  const char* ptr = infile.Line();
  if (ptr == 0) {
    mprinterr("Error: Unexpected end of XPLOR restraint file.\n");
    return 1;
  }
  while ( ptr != 0 ) {
    if (ptr[0] == 'a' && ptr[1] == 's' && ptr[2] == 's' &&
        ptr[3] == 'i' && ptr[4] == 'g' && ptr[5] == 'n'   )
    {
      // 'assign' statement
      ArgList line(ptr, " ()");
      if (line.empty()) {
        mprinterr("Error: Could not parse XPLOR 'assign' line:\n\t%s",ptr);
      } else {
        line.MarkArg(0); // Mark 'assign'
        // Get 2 Masks
        std::string maskExpression1 = GetSelection( line, resOffset_ );
        std::string maskExpression2 = GetSelection( line, resOffset_ );
        if (maskExpression1.empty() || maskExpression2.empty()) {
          mprinterr("Error: Could not get masks from line:\n\t%s", ptr);
          mprinterr("Error: Check if residue number + offset is out of bounds.\n");
        } else {
          NOE.dMask1_.SetMaskString( maskExpression1 );
          NOE.dMask2_.SetMaskString( maskExpression2 );
          // Check for noe bounds
          NOE.rexp_ = line.getNextDouble(-1.0);
          if ( NOE.rexp_ < 0.0 ) {
            // No more on this line, assume jcoupling
            ptr = infile.Line();
            line.SetList(ptr, " ()");
            // Get 2 more masks and jcoupling values 
          } else {
            // NOE
            NOE.boundh_ = NOE.rexp_ + line.getNextDouble(0.0);
            NOE.bound_ = NOE.rexp_ - line.getNextDouble(0.0);
            NOE.dist_ = 0;
            NOEs_.push_back( NOE );
          }
        }
      }
    }
    ptr = infile.Line();
  }
  infile.CloseFile();
  return 0;
}
// -----------------------------------------------------------------------------
// Action_NMRrst::Setup()
/** Determine what atoms each mask pertains to for the current parm file.
  */
Action::RetType Action_NMRrst::Setup(Topology* currentParm, Topology** parmAddress) {
  for (noeArray::iterator noe = NOEs_.begin(); noe != NOEs_.end(); ++noe) {
    if (currentParm->SetupIntegerMask( (*noe).dMask1_ )) return Action::ERR;
    if (currentParm->SetupIntegerMask( (*noe).dMask2_ )) return Action::ERR;
    //mprintf("\t%s (%i atoms) to %s (%i atoms)",Mask1_.MaskString(), Mask1_.Nselected(),
    //        Mask2_.MaskString(),Mask2_.Nselected());
    if ((*noe).dMask1_.None() || (*noe).dMask2_.None()) {
      mprintf("Warning: One or both masks for NOE '%s' have no atoms (%i and %i).\n",
              (*noe).dist_->Legend().c_str(), (*noe).dMask1_.Nselected(),
              (*noe).dMask2_.Nselected());
      (*noe).active_ = false; 
    } else
      (*noe).active_ = true;
  }
  // Set up imaging info for this parm
  Image_.SetupImaging( currentParm->BoxType() );
  if (Image_.ImagingEnabled())
    mprintf("\tImaged.\n");
  else
    mprintf("\tImaging off.\n");

  return Action::OK;  
}

// Action_NMRrst::DoAction()
Action::RetType Action_NMRrst::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  double Dist;
  Matrix_3x3 ucell, recip;
  Vec3 a1, a2;

  for (noeArray::iterator noe = NOEs_.begin(); noe != NOEs_.end(); ++noe) {
    if ( (*noe).active_ ) {
      if (useMass_) {
        a1 = currentFrame->VCenterOfMass( (*noe).dMask1_ );
        a2 = currentFrame->VCenterOfMass( (*noe).dMask2_ );
      } else {
        a1 = currentFrame->VGeometricCenter( (*noe).dMask1_ );
        a2 = currentFrame->VGeometricCenter( (*noe).dMask2_ );
      }

      switch ( Image_.ImageType() ) {
        case NONORTHO:
          currentFrame->BoxCrd().ToRecip(ucell, recip);
          Dist = DIST2_ImageNonOrtho(a1, a2, ucell, recip);
          break;
        case ORTHO:
          Dist = DIST2_ImageOrtho(a1, a2, currentFrame->BoxCrd());
          break;
        case NOIMAGE:
          Dist = DIST2_NoImage(a1, a2);
          break;
      }
      Dist = sqrt(Dist);
      (*noe).dist_->Add(frameNum, &Dist);
    }
  }
  return Action::OK;
} 
