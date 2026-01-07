#include "LJ1264_Params.h"

using namespace Cpptraj::Parm;

/** CONSTRUCTOR */
LJ1264_Params::LJ1264_Params() :
  waterModel_(UNKNOWN_WATER_MODEL),
  tunfactor_(1.0)
{}

void LJ1264_Params::set_tip3p_params() {
  c4params_.insert( NameMapPair("Li1", 27.0) );
  c4params_.insert( NameMapPair("Na1", 0.0) );
  c4params_.insert( NameMapPair("K1", 8.0) );
  c4params_.insert( NameMapPair("Rb1", 0.0) );
  c4params_.insert( NameMapPair("Cs1", 2.0) );
  c4params_.insert( NameMapPair("Tl1", 50.0) );
  c4params_.insert( NameMapPair("Cu1", 7.0) );
  c4params_.insert( NameMapPair("Ag1", 83.0) );
  c4params_.insert( NameMapPair("F-1", -27.0) );
  c4params_.insert( NameMapPair("Cl-1", -38.0) );
  c4params_.insert( NameMapPair("Br-1", -39.0) );
  c4params_.insert( NameMapPair("I-1", -45.0) );
  c4params_.insert( NameMapPair("Be2", 186.5) );
  c4params_.insert( NameMapPair("Cu2", 290.9) );
  c4params_.insert( NameMapPair("Ni2", 212.8) );
  c4params_.insert( NameMapPair("Zn2", 231.6) );
  c4params_.insert( NameMapPair("Co2", 209.7) );
  c4params_.insert( NameMapPair("Cr2", 136.8) );
  c4params_.insert( NameMapPair("Fe2", 163.0) );
  c4params_.insert( NameMapPair("Mg2", 132.9) );
  c4params_.insert( NameMapPair("V2", 195.7) );
  c4params_.insert( NameMapPair("Mn2", 146.1) );
  c4params_.insert( NameMapPair("Hg2", 288.8) );
  c4params_.insert( NameMapPair("Cd2", 185.6) );
  c4params_.insert( NameMapPair("Ca2", 87.3) );
  c4params_.insert( NameMapPair("Sn2", 187.9) );
  c4params_.insert( NameMapPair("Sr2", 82.7) );
  c4params_.insert( NameMapPair("Ba2", 71.9) );
  c4params_.insert( NameMapPair("Al3", 399.0) );
  c4params_.insert( NameMapPair("Fe3", 428.0) );
  c4params_.insert( NameMapPair("Cr3", 258.0) );
  c4params_.insert( NameMapPair("In3", 347.0) );
  c4params_.insert( NameMapPair("Tl3", 456.0) );
  c4params_.insert( NameMapPair("Y3", 216.0) );
  c4params_.insert( NameMapPair("La3", 152.0) );
  c4params_.insert( NameMapPair("Ce3", 230.0) );
  c4params_.insert( NameMapPair("Pr3", 264.0) );
  c4params_.insert( NameMapPair("Nd3", 213.0) );
  c4params_.insert( NameMapPair("Sm3", 230.0) );
  c4params_.insert( NameMapPair("Eu3", 259.0) );
  c4params_.insert( NameMapPair("Gd3", 198.0) );
  c4params_.insert( NameMapPair("Tb3", 235.0) );
  c4params_.insert( NameMapPair("Dy3", 207.0) );
  c4params_.insert( NameMapPair("Er3", 251.0) );
  c4params_.insert( NameMapPair("Tm3", 282.0) );
  c4params_.insert( NameMapPair("Lu3", 249.0) );
  c4params_.insert( NameMapPair("Hf4", 827.0) );
  c4params_.insert( NameMapPair("Zr4", 761.0) );
  c4params_.insert( NameMapPair("Ce4", 706.0) );
  c4params_.insert( NameMapPair("U4", 1034.0) );
  c4params_.insert( NameMapPair("Pu4", 828.0) );
  c4params_.insert( NameMapPair("Th4", 512.0) );
}

/** Set up C4 parameters for given water model */
void LJ1264_Params::setupForWaterModel(WaterModelType wmIn)
{
  if (wmIn == waterModel_) return;
  c4params_.clear();
  waterModel_ = wmIn;

  switch (waterModel_) {
    case UNKNOWN_WATER_MODEL : break;
    case TIP3P : set_tip3p_params(); break;
  }
}
