#include "Action_Test.h"
#include "Constants.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "Dist_Imaged.h"
#include <cmath>

// Action_Test::Help()
void Action_Test::Help() const {
  mprintf("THIS ACTION IS ONLY FOR DEBUGGING PURPOSES.\n");
}

// Action_Test::Init()
Action::RetType Action_Test::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Datasets to store distances
  //D1_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(actionArgs.GetStringNext()), "Test1");
  //D2_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(D1_->Meta().Name()), "Test2");

  outfile_ = init.DFL().AddCpptrajFile("test.cpptraj.output",
                                       "Test Cpptraj Output", DataFileList::TEXT, true);
  if (outfile_ == 0) return Action::ERR;

  mprintf("    TEST ACTION\n");
  mask1_.SetMaskString("^1");
  mprintf("\tMask1: %s\n", mask1_.MaskString());
  mask2_.SetMaskString("!^1");
  mprintf("\tMask2: %s\n", mask2_.MaskString());
  mprintf("\tOutput to %s\n", outfile_->Filename().full());
  return Action::OK;
}

// Action_Test::Setup()
Action::RetType Action_Test::Setup(ActionSetup& setup)
{
  // Determine Box info
  if (!setup.CoordInfo().TrajBox().HasBox()) {
    mprintf("Warning: Topology %s does not contain box information.\n", setup.Top().c_str());
    return Action::SKIP;
  }

  // Masks
  if (setup.Top().SetupIntegerMask( mask1_ )) return Action::ERR;
  mask1_.MaskInfo();
  if (setup.Top().SetupIntegerMask( mask2_ )) return Action::ERR;
  mask2_.MaskInfo();

  std::vector<int> molnums = setup.Top().MolnumsSelectedBy( mask2_ );
  firstAtoms_.clear();
  firstAtoms_.reserve(molnums.size());
  for (std::vector<int>::const_iterator it = molnums.begin(); it != molnums.end(); ++it)
  {
    Molecule const& mol = setup.Top().Mol( *it );
    firstAtoms_.push_back( mol.MolUnit().Front() );
  }
  return Action::OK;

}

// Action_Test::DoAction()
Action::RetType Action_Test::DoAction(int frameNum, ActionFrame& frm)
{
  Box const& box = frm.Frm().BoxCrd();
  // Determine the minimum imaged distance between molecule 1 and the first
  // atom of every other molecule.
  Vec3 mol1 = frm.Frm().VGeometricCenter( mask1_ );

  outfile_->Printf("#Frame %8i { %12.4f %12.4f %12.4f }\n", frameNum, mol1[0], mol1[1], mol1[2]);
  for (std::vector<int>::const_iterator it = firstAtoms_.begin(); it != firstAtoms_.end(); ++it)
  {
    Vec3 mol2( frm.Frm().XYZ( *it ) );
    int ixyz1[3];
    t1_.Start();
    //double d1_sq = DIST2_ImageNonOrtho(mol1.Dptr(), mol2.Dptr(), box.UnitCell(), box.FracCell());
    double d1_sq = DIST2_ImageNonOrthoRecip(box.FracCell() * mol1, box.FracCell() * mol2, -1.0, ixyz1, box.UnitCell());
    t1_.Stop();

    int ixyz2[3];
    t2_.Start();
    double d2_sq = Cpptraj::Dist2_Imaged(mol1, mol2, box.UnitCell(), box.FracCell(), ixyz2);
    t2_.Stop();
    double delta = d1_sq - d2_sq;
    if (delta < 0.0) delta = -delta;
    if (delta > Constants::SMALL ||
        ixyz1[0] != ixyz2[0] ||
        ixyz1[1] != ixyz2[1] ||
        ixyz1[2] != ixyz2[2])
      outfile_->Printf("D2 %8i %12.4f %12.4f {%2i %2i %2i} {%2i %2i %2i}\n", *it, sqrt(d1_sq), sqrt(d2_sq),
                       ixyz1[0], ixyz1[1], ixyz1[2],
                       ixyz2[0], ixyz2[1], ixyz2[2]);

    int ixyz3[3];
    t3_.Start();
    double d3_sq = Cpptraj::Dist2_Imaged_Cart(mol1, mol2, box.UnitCell(), box.FracCell(), ixyz3);
    t3_.Stop();
    delta = d1_sq - d3_sq;
    if (delta < 0.0) delta = -delta;
    if (delta > Constants::SMALL ||
        ixyz1[0] != ixyz3[0] ||
        ixyz1[1] != ixyz3[1] ||
        ixyz1[2] != ixyz3[2])
      outfile_->Printf("D3 %8i %12.4f %12.4f {%2i %2i %2i} {%2i %2i %2i}\n", *it, sqrt(d1_sq), sqrt(d3_sq),
                       ixyz1[0], ixyz1[1], ixyz1[2],
                       ixyz3[0], ixyz3[1], ixyz3[2]);
  }
  return Action::OK; 
}

void Action_Test::Print() {
  mprintf("    TEST ACTION.\n");
  t1_.WriteTiming(1, "Original dist :");
  t2_.WriteTiming(1, "New dist      :");
  t3_.WriteTiming(1, "New dist A    :");
}
