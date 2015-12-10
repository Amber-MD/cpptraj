#include "Action_MultiVector.h"
#include "CpptrajStdio.h"

Action_MultiVector::Action_MultiVector() :
  debug_(0),
  outfile_(0),
  masterDSL_(0),
  ired_(false)
{}

void Action_MultiVector::Help() const {
  mprintf("\t[<name>] [resrange <range>] name1 <name1> name2 <name2> [out <filename>]\n"
          "\t[ired]\n"
          "  Calculate vectors between named atoms for residues in given <range>.\n");
}

inline static int SetName(NameType& name, std::string const& expr, const char* str) {
  if (expr.empty()) {
    mprinterr("Error: %s not specified.\n", str);
    return 1;
  }
  name = expr;
  return 0;
}

Action::RetType Action_MultiVector::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  // Get keywords
  outfile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs);
  std::string resrange_arg = actionArgs.GetStringKey("resrange");
  if (!resrange_arg.empty())
    if (resRange_.SetRange( resrange_arg )) return Action::ERR;
  ired_ = actionArgs.hasKey("ired");
  // Get atom names
  if (SetName(name1_, actionArgs.GetStringKey("name1"), "name1")) return Action::ERR;
  if (SetName(name2_, actionArgs.GetStringKey("name2"), "name2")) return Action::ERR;
  // Setup DataSet(s) name
  dsetname_ = actionArgs.GetStringNext();

  mprintf("    MULTIVECTOR: Calculating");
  if (ired_) mprintf(" IRED");
  if (!resRange_.Empty())
    mprintf(" vectors for residues in range %s\n", resRange_.RangeArg());
  else
    mprintf(" vectors for all solute residues.\n");
  mprintf("\tName1='%s' (origin)  Name2='%s'\n", *name1_, *name2_);
  if (!dsetname_.empty())
    mprintf("\tDataSet name: %s\n", dsetname_.c_str());
  if (outfile_ != 0) mprintf("\tOutput to %s\n", outfile_->DataFilename().base());
  init.DSL().SetDataSetsPending(true);
  masterDSL_ = init.DslPtr();
  return Action::OK;
}

// Action_MultiVector::Setup();
Action::RetType Action_MultiVector::Setup(ActionSetup& setup) {
  Range actualRange;
  // If range is empty (i.e. no resrange arg given) look through all 
  // solute residues.
  if (resRange_.Empty())
    actualRange = setup.Top().SoluteResidues();
  else {
    // If user range specified, create new range shifted by -1 since internal
    // resnums start from 0.
    actualRange = resRange_;
    actualRange.ShiftBy(-1);
  }
  // Exit if no residues specified
  if (actualRange.Empty()) {
    mprintf("Warning: No residues specified for %s\n",setup.Top().c_str());
    return Action::SKIP;
  }
  // Set default DataSet name if not specified.
  if (dsetname_.empty())
    dsetname_ = masterDSL_->GenerateDefaultName( "MVEC" );
  // Search for specified atom names in each residue in the range
  CrdIdx1_.clear();
  CrdIdx2_.clear();
  for (Range::const_iterator res = actualRange.begin(); res != actualRange.end(); ++res)
  {
    int atom1 = setup.Top().FindAtomInResidue( *res, name1_ );
    int atom2 = setup.Top().FindAtomInResidue( *res, name2_ );
    if (atom1 != -1 && atom2 != -1) {
      MetaData md(dsetname_, atom1+1);
      md.SetScalarMode( MetaData::M_VECTOR );
      if (ired_) md.SetScalarType( MetaData::IREDVEC );
      DataSet_Vector* ds = (DataSet_Vector*)masterDSL_->CheckForSet( md );
      if (ds == 0) {
        // Create DataSet
        ds = (DataSet_Vector*)masterDSL_->AddSet( DataSet::VECTOR, md );
        if (ds == 0) return Action::ERR;
        ds->SetLegend( "v" + setup.Top().AtomMaskName(atom1) + "->" +
                             setup.Top().AtomMaskName(atom2) );
        if (outfile_ != 0) outfile_->AddDataSet( ds );
      }
      data_.push_back( ds );
      CrdIdx1_.push_back( atom1 * 3 ); // Pre calc coordinate index
      CrdIdx2_.push_back( atom2 * 3 );
    } else if ((atom1==-1) != (atom2==-1)) {
      if (atom1==-1)
        mprintf("Warning: '%s' not found but '%s' found for residue %i.\n", 
                *name1_, *name2_, *res + 1);
      else // atom2==-1
        mprintf("Warning: '%s' not found but '%s' found for residue %i.\n",
                *name2_, *name1_, *res + 1);
    }
  }
  mprintf("\tSelected %zu vectors.\n", CrdIdx1_.size());
  for (std::vector<DataSet_Vector*>::const_iterator it = data_.begin();
                                                    it != data_.end(); ++it)
    mprintf("\t  %s\n", (*it)->legend());

  return Action::OK;
}

// Action_MultiVector::DoAction()
Action::RetType Action_MultiVector::DoAction(int frameNum, ActionFrame& frm) {
  for (unsigned int nv = 0; nv < CrdIdx1_.size(); ++nv) {
    Vec3 CXYZ( frm.Frm().CRD( CrdIdx1_[nv] ) );
    Vec3 VXYZ( frm.Frm().CRD( CrdIdx2_[nv] ) );
    VXYZ -= CXYZ;
    data_[nv]->AddVxyz(VXYZ, CXYZ);
  }

  return Action::OK;
}
