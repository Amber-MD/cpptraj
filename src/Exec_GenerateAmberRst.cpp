#include "Exec_GenerateAmberRst.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "TorsionRoutines.h"
#include "Constants.h"

void Exec_GenerateAmberRst::Help() const {
  mprintf("\t<mask1> <mask2> [<mask3>] [<mask4>]\n"
          "\tr1 <r1> r2 <r2> r3 <r3> r4 <r4> rk2 <rk2> rk3 <rk3>\n"
          "\t{%s}\n"
          "\t[{%s} [offset <off>] [width <width>]]\n"
          "\t[out <outfile>]\n"
          "  Generate Amber-format restraint from 2 or more mask expressions.\n",
          DataSetList::TopArgs, DataSetList::RefArgs);
}

/// Generate amber restraints from given masks.
Exec::RetType Exec_GenerateAmberRst::Execute(CpptrajState& State, ArgList& argIn)
{
  // Get parm
  Topology* parm = State.DSL().GetTopology( argIn );
  if (parm == 0) {
    mprinterr("Error: No parm files loaded.\n");
    return CpptrajState::ERR;
  }
  // Get optional reference coords
  ReferenceFrame RefCrd = State.DSL().GetReferenceFrame(argIn);
  // Get arguments
  if (argIn.hasKey("overwrite"))
    mprintf("Warning: 'overwrite' keyword no longer necessary and is deprecated.\n");
  double r1 = argIn.getKeyDouble("r1", 0.0);
  double r2 = argIn.getKeyDouble("r2", 0.0);
  double r3 = argIn.getKeyDouble("r3", 0.0);
  double r4 = argIn.getKeyDouble("r4", 0.0);
  double rk2 = argIn.getKeyDouble("rk2", 0.0);
  double rk3 = argIn.getKeyDouble("rk3", 0.0);
  // crddist will be !RefCrd.empty()
  double offset = argIn.getKeyDouble("offset", 0.0);
  double width = argIn.getKeyDouble("width", 0.5);
  CpptrajFile* outfile = State.DFL().AddCpptrajFile(argIn.GetStringKey("out"), "Amber Rst",
                                                    DataFileList::TEXT, true);
  if (outfile == 0) {
    mprinterr("Error: Could not open output file.\n");
    return CpptrajState::ERR;
  }
  // TODO: else if (strcmp(pch,"nobb")==0) nobb=1;
  // Assume everything else is a mask
  std::vector<AtomMask> rstMasks;
  std::string maskExpr = argIn.GetMaskNext();
  while (!maskExpr.empty()) {
    if ( rstMasks.size() >= 4 )
      mprintf("Warning: 4 masks already defined. Skipping '%s'\n", maskExpr.c_str());
    else {
      AtomMask tmpmask( maskExpr );
      int maskerr = 0;
      if (!RefCrd.empty())
        maskerr = parm->SetupIntegerMask( tmpmask, RefCrd.Coord() );
      else
        maskerr = parm->SetupIntegerMask( tmpmask );
      if ( maskerr != 0 ) {
        mprinterr("Error: Could not set up mask '%s'\n", tmpmask.MaskString());
        return CpptrajState::ERR;
      }
      tmpmask.MaskInfo();
      if ( tmpmask.None() ) {
        mprinterr("Error: '%s' corresponds to no atoms.\n", tmpmask.MaskString());
        return CpptrajState::ERR;
      }
      rstMasks.push_back( tmpmask );
    }
    maskExpr = argIn.GetMaskNext();
  }
  if (State.Debug() > 0) mprintf("\tDefined %zu rst masks.\n", rstMasks.size());
  if (rstMasks.size() < 2) {
    mprinterr("Error: Must specify at least 2 masks for restraint.\n");
    return CpptrajState::ERR;
  }
  // TODO: Remove backbone atoms for 'nobb'?
  // If a reference frame was specified and distance restraint, use center of
  // mass distance/angle/torsion between masks as r2.
  if ( !RefCrd.empty() ) {
    if ( RefCrd.Parm().Pindex() != parm->Pindex() )
      mprintf("Warning: Reference topology does not match specified topology.\n");
    Vec3 a1 = RefCrd.Coord().VCenterOfMass( rstMasks[0] );
    Vec3 a2 = RefCrd.Coord().VCenterOfMass( rstMasks[1] );
    if (rstMasks.size() == 2)
      r2 = DIST_NoImage( a1, a2 );
    else if (rstMasks.size() == 3) {
      Vec3 a3 = RefCrd.Coord().VCenterOfMass( rstMasks[2] );
      r2 = CalcAngle(a1.Dptr(), a2.Dptr(), a3.Dptr()) * Constants::RADDEG;
    } else if (rstMasks.size() == 4) {
      Vec3 a3 = RefCrd.Coord().VCenterOfMass( rstMasks[2] );
      Vec3 a4 = RefCrd.Coord().VCenterOfMass( rstMasks[3] );
      r2 = Torsion(a1.Dptr(), a2.Dptr(), a3.Dptr(), a4.Dptr()) * Constants::RADDEG;
    }
    r2 += offset;
    r3 = r2;
    r1 = r2 - width;
    r4 = r3 + width;
    mprintf("\tCoM value from ref will be used, r1=%f, r2=%f, r3=%f, r4=%f\n", r1,r2,r3,r4);
  }
  // Print restraint header 
  outfile->Printf(" &rst iat=");
  for (std::vector<AtomMask>::const_iterator M = rstMasks.begin();
                                             M != rstMasks.end(); ++M)
  {
    if ((*M).Nselected() == 1)
      outfile->Printf("%i,", (*M)[0] + 1);
    else
      outfile->Printf("-1,");
  }
  outfile->Printf("0\n");
  // Print Restraint boundaries and constants
  outfile->Printf("   r1=%f, r2=%f, r3=%f, r4=%f, rk2=%f, rk3=%f,\n",
                 r1, r2, r3, r4, rk2, rk3);
  // Print out atom groups if necessary
  unsigned int group = 1;
  for (std::vector<AtomMask>::const_iterator M = rstMasks.begin();
                                             M != rstMasks.end(); ++M, group++)
  {
    if ((*M).Nselected() > 1) {
      outfile->Printf("   ");
      unsigned int j = 1;
      for (AtomMask::const_iterator atom = (*M).begin();
                                    atom != (*M).end(); ++atom, j++)
        outfile->Printf("IGR%u(%u)=%i,", group, j, (*atom) + 1);
      outfile->Printf("\n");
    }
  }
  // Finish restraint
  outfile->Printf("   nstep1=0, nstep2=0,\n &end\n");
  return CpptrajState::OK;
}
