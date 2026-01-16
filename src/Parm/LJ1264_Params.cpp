#include "LJ1264_Params.h"
#include "../ArgList.h"
#include "../BufferedLine.h"
#include "../Constants.h" // SMALL
#include "../CpptrajStdio.h"
#include "../StringRoutines.h" // integerToString
#include "../Topology.h"
#include <cstdlib> // atof
#include <cctype> // tolower
#include <cmath>

using namespace Cpptraj::Parm;

const double LJ1264_Params::DEFAULT_WATER_POL_ = 1.444; ///< Polarizability of water

/** CONSTRUCTOR */
LJ1264_Params::LJ1264_Params() :
  waterModel_(UNKNOWN_WATER_MODEL),
  tunfactor_(1.0),
  WATER_POL_(DEFAULT_WATER_POL_),
  debug_(0)
{}

// ----- WATER MODEL LJC PARAMETERS GO HERE ------------------------------------
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

void LJ1264_Params::set_tip4pew_params() {
  c4params_.insert( NameMapPair("Li1", 36.0) );
  c4params_.insert( NameMapPair("Na1", 9.0) );
  c4params_.insert( NameMapPair("K1", 24.0) );
  c4params_.insert( NameMapPair("Rb1", 13.0) );
  c4params_.insert( NameMapPair("Cs1", 16.0) );
  c4params_.insert( NameMapPair("Tl1", 65.0) );
  c4params_.insert( NameMapPair("Cu1", 21.0) );
  c4params_.insert( NameMapPair("Ag1", 94.0) );
  c4params_.insert( NameMapPair("F-1", -67.0) );
  c4params_.insert( NameMapPair("Cl-1", -66.0) );
  c4params_.insert( NameMapPair("Br-1", -68.0) );
  c4params_.insert( NameMapPair("I-1", -62.0) );
  c4params_.insert( NameMapPair("Be2", 228.5) );
  c4params_.insert( NameMapPair("Cu2", 339.2) );
  c4params_.insert( NameMapPair("Ni2", 259.2) );
  c4params_.insert( NameMapPair("Zn2", 272.3) );
  c4params_.insert( NameMapPair("Co2", 252.8) );
  c4params_.insert( NameMapPair("Cr2", 177.4) );
  c4params_.insert( NameMapPair("Fe2", 201.1) );
  c4params_.insert( NameMapPair("Mg2", 180.5) );
  c4params_.insert( NameMapPair("V2", 244.8) );
  c4params_.insert( NameMapPair("Mn2", 192.3) );
  c4params_.insert( NameMapPair("Hg2", 335.2) );
  c4params_.insert( NameMapPair("Cd2", 233.7) );
  c4params_.insert( NameMapPair("Ca2", 128.0) );
  c4params_.insert( NameMapPair("Sn2", 231.4) );
  c4params_.insert( NameMapPair("Sr2", 118.9) );
  c4params_.insert( NameMapPair("Ba2", 112.5) );
  c4params_.insert( NameMapPair("Al3", 488.0) );
  c4params_.insert( NameMapPair("Fe3", 519.0) );
  c4params_.insert( NameMapPair("Cr3", 322.0) );
  c4params_.insert( NameMapPair("In3", 425.0) );
  c4params_.insert( NameMapPair("Tl3", 535.0) );
  c4params_.insert( NameMapPair("Y3", 294.0) );
  c4params_.insert( NameMapPair("La3", 243.0) );
  c4params_.insert( NameMapPair("Ce3", 315.0) );
  c4params_.insert( NameMapPair("Pr3", 348.0) );
  c4params_.insert( NameMapPair("Nd3", 297.0) );
  c4params_.insert( NameMapPair("Sm3", 314.0) );
  c4params_.insert( NameMapPair("Eu3", 345.0) );
  c4params_.insert( NameMapPair("Gd3", 280.0) );
  c4params_.insert( NameMapPair("Tb3", 313.0) );
  c4params_.insert( NameMapPair("Dy3", 298.0) );
  c4params_.insert( NameMapPair("Er3", 328.0) );
  c4params_.insert( NameMapPair("Tm3", 356.0) );
  c4params_.insert( NameMapPair("Lu3", 331.0) );
  c4params_.insert( NameMapPair("Hf4", 956.0) );
  c4params_.insert( NameMapPair("Zr4", 895.0) );
  c4params_.insert( NameMapPair("Ce4", 835.0) );
  c4params_.insert( NameMapPair("U4", 1183.0) );
  c4params_.insert( NameMapPair("Pu4", 972.0) );
  c4params_.insert( NameMapPair("Th4", 625.0) );
}

void LJ1264_Params::set_spce_params() {
  c4params_.insert( NameMapPair("Li1", 33.0) );
  c4params_.insert( NameMapPair("Na1", 6.0) );
  c4params_.insert( NameMapPair("K1", 19.0) );
  c4params_.insert( NameMapPair("Rb1", 7.0) );
  c4params_.insert( NameMapPair("Cs1", 12.0) );
  c4params_.insert( NameMapPair("Tl1", 61.0) );
  c4params_.insert( NameMapPair("Cu1", 9.0) );
  c4params_.insert( NameMapPair("Ag1", 92.0) );
  c4params_.insert( NameMapPair("F-1", -53.0) );
  c4params_.insert( NameMapPair("Cl-1", -55.0) );
  c4params_.insert( NameMapPair("Br-1", -51.0) );
  c4params_.insert( NameMapPair("I-1", -51.0) );
  c4params_.insert( NameMapPair("Be2", 188.1) );
  c4params_.insert( NameMapPair("Cu2", 304.4) );
  c4params_.insert( NameMapPair("Ni2", 205.2) );
  c4params_.insert( NameMapPair("Zn2", 231.2) );
  c4params_.insert( NameMapPair("Co2", 209.2) );
  c4params_.insert( NameMapPair("Cr2", 131.2) );
  c4params_.insert( NameMapPair("Fe2", 155.4) );
  c4params_.insert( NameMapPair("Mg2", 122.2) );
  c4params_.insert( NameMapPair("V2", 206.6) );
  c4params_.insert( NameMapPair("Mn2", 154.9) );
  c4params_.insert( NameMapPair("Hg2", 300.2) );
  c4params_.insert( NameMapPair("Cd2", 198.8) );
  c4params_.insert( NameMapPair("Ca2", 89.0) );
  c4params_.insert( NameMapPair("Sn2", 201.1) );
  c4params_.insert( NameMapPair("Sr2", 96.3) );
  c4params_.insert( NameMapPair("Ba2", 85.8) );
  c4params_.insert( NameMapPair("Al3", 406.0) );
  c4params_.insert( NameMapPair("Fe3", 442.0) );
  c4params_.insert( NameMapPair("Cr3", 254.0) );
  c4params_.insert( NameMapPair("In3", 349.0) );
  c4params_.insert( NameMapPair("Tl3", 455.0) );
  c4params_.insert( NameMapPair("Y3", 209.0) );
  c4params_.insert( NameMapPair("La3", 165.0) );
  c4params_.insert( NameMapPair("Ce3", 242.0) );
  c4params_.insert( NameMapPair("Pr3", 272.0) );
  c4params_.insert( NameMapPair("Nd3", 235.0) );
  c4params_.insert( NameMapPair("Sm3", 224.0) );
  c4params_.insert( NameMapPair("Eu3", 273.0) );
  c4params_.insert( NameMapPair("Gd3", 186.0) );
  c4params_.insert( NameMapPair("Tb3", 227.0) );
  c4params_.insert( NameMapPair("Dy3", 206.0) );
  c4params_.insert( NameMapPair("Er3", 247.0) );
  c4params_.insert( NameMapPair("Tm3", 262.0) );
  c4params_.insert( NameMapPair("Lu3", 247.0) );
  c4params_.insert( NameMapPair("Hf4", 810.0) );
  c4params_.insert( NameMapPair("Zr4", 760.0) );
  c4params_.insert( NameMapPair("Ce4", 694.0) );
  c4params_.insert( NameMapPair("U4", 1043.0) );
  c4params_.insert( NameMapPair("Pu4", 828.0) );
  c4params_.insert( NameMapPair("Th4", 513.0) );
}

void LJ1264_Params::set_opc3_params() {
  c4params_.insert( NameMapPair("Li1", 29.0) );
  c4params_.insert( NameMapPair("Na1", 2.0) );
  c4params_.insert( NameMapPair("K1", 16.0) );
  c4params_.insert( NameMapPair("Rb1", 8.0) );
  c4params_.insert( NameMapPair("Cs1", 6.0) );
  c4params_.insert( NameMapPair("Tl1", 63.0) );
  c4params_.insert( NameMapPair("Cu1", 12.0) );
  c4params_.insert( NameMapPair("Ag1", 83.0) );
  c4params_.insert( NameMapPair("F-1", -40.0) );
  c4params_.insert( NameMapPair("Cl-1", -47.0) );
  c4params_.insert( NameMapPair("Br-1", -43.0) );
  c4params_.insert( NameMapPair("I-1", -45.0) );
  c4params_.insert( NameMapPair("Be2", 186.0) );
  c4params_.insert( NameMapPair("Cu2", 269.0) );
  c4params_.insert( NameMapPair("Ni2", 207.0) );
  c4params_.insert( NameMapPair("Zn2", 199.0) );
  c4params_.insert( NameMapPair("Co2", 182.0) );
  c4params_.insert( NameMapPair("Cr2", 109.0) );
  c4params_.insert( NameMapPair("Fe2", 131.0) );
  c4params_.insert( NameMapPair("Mg2", 117.0) );
  c4params_.insert( NameMapPair("V2", 201.0) );
  c4params_.insert( NameMapPair("Mn2", 137.0) );
  c4params_.insert( NameMapPair("Hg2", 276.0) );
  c4params_.insert( NameMapPair("Cd2", 200.0) );
  c4params_.insert( NameMapPair("Ca2", 76.0) );
  c4params_.insert( NameMapPair("Sn2", 188.0) );
  c4params_.insert( NameMapPair("Sr2", 85.0) );
  c4params_.insert( NameMapPair("Ba2", 77.0) );
  c4params_.insert( NameMapPair("Al3", 363.0) );
  c4params_.insert( NameMapPair("Fe3", 429.0) );
  c4params_.insert( NameMapPair("Cr3", 209.0) );
  c4params_.insert( NameMapPair("In3", 330.0) );
  c4params_.insert( NameMapPair("Tl3", 437.0) );
  c4params_.insert( NameMapPair("Y3", 192.0) );
  c4params_.insert( NameMapPair("La3", 131.0) );
  c4params_.insert( NameMapPair("Ce3", 215.0) );
  c4params_.insert( NameMapPair("Pr3", 255.0) );
  c4params_.insert( NameMapPair("Nd3", 184.0) );
  c4params_.insert( NameMapPair("Sm3", 188.0) );
  c4params_.insert( NameMapPair("Eu3", 233.0) );
  c4params_.insert( NameMapPair("Gd3", 164.0) );
  c4params_.insert( NameMapPair("Tb3", 199.0) );
  c4params_.insert( NameMapPair("Dy3", 183.0) );
  c4params_.insert( NameMapPair("Er3", 228.0) );
  c4params_.insert( NameMapPair("Tm3", 246.0) );
  c4params_.insert( NameMapPair("Lu3", 222.0) );
  c4params_.insert( NameMapPair("Hf4", 718.0) );
  c4params_.insert( NameMapPair("Zr4", 707.0) );
  c4params_.insert( NameMapPair("Ce4", 653.0) );
  c4params_.insert( NameMapPair("U4", 980.0) );
  c4params_.insert( NameMapPair("Pu4", 817.0) );
  c4params_.insert( NameMapPair("Th4", 452.0) );
}

void LJ1264_Params::set_opc_params() {
  c4params_.insert( NameMapPair("Li1", 29.0) );
  c4params_.insert( NameMapPair("Na1", 0.0) );
  c4params_.insert( NameMapPair("K1", 20.0) );
  c4params_.insert( NameMapPair("Rb1", 6.0) );
  c4params_.insert( NameMapPair("Cs1", 13.0) );
  c4params_.insert( NameMapPair("Tl1", 60.0) );
  c4params_.insert( NameMapPair("Cu1", 16.0) );
  c4params_.insert( NameMapPair("Ag1", 83.0) );
  c4params_.insert( NameMapPair("F-1", -67.0) );
  c4params_.insert( NameMapPair("Cl-1", -69.0) );
  c4params_.insert( NameMapPair("Br-1", -60.0) );
  c4params_.insert( NameMapPair("I-1", -60.0) );
  c4params_.insert( NameMapPair("Be2", 214.0) );
  c4params_.insert( NameMapPair("Cu2", 291.0) );
  c4params_.insert( NameMapPair("Ni2", 212.0) );
  c4params_.insert( NameMapPair("Zn2", 225.0) );
  c4params_.insert( NameMapPair("Co2", 204.0) );
  c4params_.insert( NameMapPair("Cr2", 132.0) );
  c4params_.insert( NameMapPair("Fe2", 154.0) );
  c4params_.insert( NameMapPair("Mg2", 127.0) );
  c4params_.insert( NameMapPair("V2", 239.0) );
  c4params_.insert( NameMapPair("Mn2", 175.0) );
  c4params_.insert( NameMapPair("Hg2", 289.0) );
  c4params_.insert( NameMapPair("Cd2", 219.0) );
  c4params_.insert( NameMapPair("Ca2", 86.0) );
  c4params_.insert( NameMapPair("Sn2", 199.0) );
  c4params_.insert( NameMapPair("Sr2", 87.0) );
  c4params_.insert( NameMapPair("Ba2", 78.0) );
  c4params_.insert( NameMapPair("Al3", 399.0) );
  c4params_.insert( NameMapPair("Fe3", 531.0) );
  c4params_.insert( NameMapPair("Cr3", 243.0) );
  c4params_.insert( NameMapPair("In3", 413.0) );
  c4params_.insert( NameMapPair("Tl3", 479.0) );
  c4params_.insert( NameMapPair("Y3", 260.0) );
  c4params_.insert( NameMapPair("La3", 165.0) );
  c4params_.insert( NameMapPair("Ce3", 289.0) );
  c4params_.insert( NameMapPair("Pr3", 311.0) );
  c4params_.insert( NameMapPair("Nd3", 243.0) );
  c4params_.insert( NameMapPair("Sm3", 236.0) );
  c4params_.insert( NameMapPair("Eu3", 279.0) );
  c4params_.insert( NameMapPair("Gd3", 222.0) );
  c4params_.insert( NameMapPair("Tb3", 256.0) );
  c4params_.insert( NameMapPair("Dy3", 243.0) );
  c4params_.insert( NameMapPair("Er3", 298.0) );
  c4params_.insert( NameMapPair("Tm3", 314.0) );
  c4params_.insert( NameMapPair("Lu3", 289.0) );
  c4params_.insert( NameMapPair("Hf4", 847.0) );
  c4params_.insert( NameMapPair("Zr4", 804.0) );
  c4params_.insert( NameMapPair("Ce4", 789.0) );
  c4params_.insert( NameMapPair("U4", 1123.0) );
  c4params_.insert( NameMapPair("Pu4", 941.0) );
  c4params_.insert( NameMapPair("Th4", 598.0) );
}

void LJ1264_Params::set_fb3_params() {
  c4params_.insert( NameMapPair("Li1", 30.0) );
  c4params_.insert( NameMapPair("Na1", 2.0) );
  c4params_.insert( NameMapPair("K1", 15.0) );
  c4params_.insert( NameMapPair("Rb1", 7.0) );
  c4params_.insert( NameMapPair("Cs1", 17.0) );
  c4params_.insert( NameMapPair("Tl1", 65.0) );
  c4params_.insert( NameMapPair("Cu1", 17.0) );
  c4params_.insert( NameMapPair("Ag1", 85.0) );
  c4params_.insert( NameMapPair("F-1", -45.0) );
  c4params_.insert( NameMapPair("Cl-1", -49.0) );
  c4params_.insert( NameMapPair("Br-1", -40.0) );
  c4params_.insert( NameMapPair("I-1", -52.0) );
  c4params_.insert( NameMapPair("Be2", 193.0) );
  c4params_.insert( NameMapPair("Cu2", 279.0) );
  c4params_.insert( NameMapPair("Ni2", 223.0) );
  c4params_.insert( NameMapPair("Zn2", 217.0) );
  c4params_.insert( NameMapPair("Co2", 192.0) );
  c4params_.insert( NameMapPair("Cr2", 138.0) );
  c4params_.insert( NameMapPair("Fe2", 157.0) );
  c4params_.insert( NameMapPair("Mg2", 128.0) );
  c4params_.insert( NameMapPair("V2", 212.0) );
  c4params_.insert( NameMapPair("Mn2", 149.0) );
  c4params_.insert( NameMapPair("Hg2", 289.0) );
  c4params_.insert( NameMapPair("Cd2", 201.0) );
  c4params_.insert( NameMapPair("Ca2", 92.0) );
  c4params_.insert( NameMapPair("Sn2", 205.0) );
  c4params_.insert( NameMapPair("Sr2", 91.0) );
  c4params_.insert( NameMapPair("Ba2", 78.0) );
  c4params_.insert( NameMapPair("Al3", 387.0) );
  c4params_.insert( NameMapPair("Fe3", 446.0) );
  c4params_.insert( NameMapPair("Cr3", 232.0) );
  c4params_.insert( NameMapPair("In3", 343.0) );
  c4params_.insert( NameMapPair("Tl3", 464.0) );
  c4params_.insert( NameMapPair("Y3", 218.0) );
  c4params_.insert( NameMapPair("La3", 155.0) );
  c4params_.insert( NameMapPair("Ce3", 251.0) );
  c4params_.insert( NameMapPair("Pr3", 291.0) );
  c4params_.insert( NameMapPair("Nd3", 211.0) );
  c4params_.insert( NameMapPair("Sm3", 218.0) );
  c4params_.insert( NameMapPair("Eu3", 261.0) );
  c4params_.insert( NameMapPair("Gd3", 191.0) );
  c4params_.insert( NameMapPair("Tb3", 226.0) );
  c4params_.insert( NameMapPair("Dy3", 219.0) );
  c4params_.insert( NameMapPair("Er3", 256.0) );
  c4params_.insert( NameMapPair("Tm3", 285.0) );
  c4params_.insert( NameMapPair("Lu3", 250.0) );
  c4params_.insert( NameMapPair("Hf4", 787.0) );
  c4params_.insert( NameMapPair("Zr4", 769.0) );
  c4params_.insert( NameMapPair("Ce4", 701.0) );
  c4params_.insert( NameMapPair("U4", 1044.0) );
  c4params_.insert( NameMapPair("Pu4", 817.0) );
  c4params_.insert( NameMapPair("Th4", 507.0) );
}

void LJ1264_Params::set_fb4_params() {
  c4params_.insert( NameMapPair("Li1", 33.0) );
  c4params_.insert( NameMapPair("Na1", 8.0) );
  c4params_.insert( NameMapPair("K1", 25.0) );
  c4params_.insert( NameMapPair("Rb1", 9.0) );
  c4params_.insert( NameMapPair("Cs1", 13.0) );
  c4params_.insert( NameMapPair("Tl1", 68.0) );
  c4params_.insert( NameMapPair("Cu1", 25.0) );
  c4params_.insert( NameMapPair("Ag1", 90.0) );
  c4params_.insert( NameMapPair("F-1", -57.0) );
  c4params_.insert( NameMapPair("Cl-1", -55.0) );
  c4params_.insert( NameMapPair("Br-1", -51.0) );
  c4params_.insert( NameMapPair("I-1", -53.0) );
  c4params_.insert( NameMapPair("Be2", 227.0) );
  c4params_.insert( NameMapPair("Cu2", 313.0) );
  c4params_.insert( NameMapPair("Ni2", 218.0) );
  c4params_.insert( NameMapPair("Zn2", 239.0) );
  c4params_.insert( NameMapPair("Co2", 206.0) );
  c4params_.insert( NameMapPair("Cr2", 159.0) );
  c4params_.insert( NameMapPair("Fe2", 187.0) );
  c4params_.insert( NameMapPair("Mg2", 133.0) );
  c4params_.insert( NameMapPair("V2", 234.0) );
  c4params_.insert( NameMapPair("Mn2", 181.0) );
  c4params_.insert( NameMapPair("Hg2", 331.0) );
  c4params_.insert( NameMapPair("Cd2", 227.0) );
  c4params_.insert( NameMapPair("Ca2", 109.0) );
  c4params_.insert( NameMapPair("Sn2", 215.0) );
  c4params_.insert( NameMapPair("Sr2", 103.0) );
  c4params_.insert( NameMapPair("Ba2", 95.0) );
  c4params_.insert( NameMapPair("Al3", 427.0) );
  c4params_.insert( NameMapPair("Fe3", 502.0) );
  c4params_.insert( NameMapPair("Cr3", 286.0) );
  c4params_.insert( NameMapPair("In3", 403.0) );
  c4params_.insert( NameMapPair("Tl3", 514.0) );
  c4params_.insert( NameMapPair("Y3", 268.0) );
  c4params_.insert( NameMapPair("La3", 211.0) );
  c4params_.insert( NameMapPair("Ce3", 294.0) );
  c4params_.insert( NameMapPair("Pr3", 326.0) );
  c4params_.insert( NameMapPair("Nd3", 256.0) );
  c4params_.insert( NameMapPair("Sm3", 257.0) );
  c4params_.insert( NameMapPair("Eu3", 302.0) );
  c4params_.insert( NameMapPair("Gd3", 238.0) );
  c4params_.insert( NameMapPair("Tb3", 270.0) );
  c4params_.insert( NameMapPair("Dy3", 253.0) );
  c4params_.insert( NameMapPair("Er3", 304.0) );
  c4params_.insert( NameMapPair("Tm3", 320.0) );
  c4params_.insert( NameMapPair("Lu3", 303.0) );
  c4params_.insert( NameMapPair("Hf4", 837.0) );
  c4params_.insert( NameMapPair("Zr4", 845.0) );
  c4params_.insert( NameMapPair("Ce4", 771.0) );
  c4params_.insert( NameMapPair("U4", 1140.0) );
  c4params_.insert( NameMapPair("Pu4", 919.0) );
  c4params_.insert( NameMapPair("Th4", 601.0) );
}

// -----------------------------------------------------------------------------

/** Set up default C4 parameters for given water model */
void LJ1264_Params::setupForWaterModel(WaterModelType wmIn)
{
  if (wmIn == waterModel_) return;
  c4params_.clear();
  waterModel_ = wmIn;

  switch (waterModel_) {
    case UNKNOWN_WATER_MODEL : return;
    case TIP3P   : set_tip3p_params(); break;
    case TIP4PEW : set_tip4pew_params(); break;
    case SPCE    : set_spce_params(); break;
    case OPC3    : set_opc3_params(); break;
    case OPC     : set_opc_params(); break;
    case FB3     : set_fb3_params(); break;
    case FB4     : set_fb4_params(); break;
  }

  mprintf("\t  Using default C4 parameters for %s\n", WaterModelStr(waterModel_));
}

/** Read a 2 column file with format <name> <value> */
int LJ1264_Params::read_2col_file(std::string const& infileName, NameMapType& mapIn)
{
  mapIn.clear();
  BufferedLine infile;
  if (infile.OpenFileRead(infileName)) {
    mprinterr("Error: Could not open '%s'\n", infileName.c_str());
    return 1;
  }
  const char* ptr = infile.Line();
  while (ptr != 0) {
    ArgList line(ptr, " \t\"'");
    if (line.Nargs() > 0 && line[0][0] != '#')
    {
      if (line.Nargs() < 2) {
        mprinterr("Error: Expected a file with at least 2 columns. Line %i has %i: %s\n",
                  infile.LineNumber(), line.Nargs(), ptr);
        return 1;
      }
      double col2 = atof( line[1].c_str() );
      // TODO report doubles
      mapIn.insert( NameMapPair(line[0], col2) );
    }
    ptr = infile.Line();
  }
  infile.CloseFile();
  return 0;
}

/** Read polarizabilities. */
int LJ1264_Params::read_pol(std::string const& polfile)
{
  mprintf("\t  Reading atomic polarizabilities from '%s'\n", polfile.c_str());
  return read_2col_file(polfile, pol_);
}

/** Read C4 params */
int LJ1264_Params::read_c4(std::string const& c4file)
{
  mprintf("\t  Reading C4 parameters from '%s'\n", c4file.c_str());
  waterModel_ = UNKNOWN_WATER_MODEL;
  return read_2col_file(c4file, c4params_);
}

/** Initialize */
int LJ1264_Params::Init_LJ1264(std::string const& maskIn, std::string const& c4fileIn, WaterModelType wmIn,
                              std::string const& polfileIn, double tunfactorIn, int debugIn)
{
  mprintf("\tEnabling LJ 12-6-4 parameters.\n");
  debug_ = debugIn;
  tunfactor_ = tunfactorIn;
  mprintf("\t  LJ 12-6-4 tuning factor: %g\n", tunfactor_);

  if (maskIn.empty()) {
    mprintf("Warning: No LJ 12-6-4 mask specified; using default.\n");
    mask_.SetMaskString(":ZN");
  } else {
    if (mask_.SetMaskString(maskIn)) {
      mprinterr("Error: Could not set mask string '%s'\n", maskIn.c_str());
      return 1;
    }
  }
  mprintf("\t  LJ 12-6-4 mask: %s\n", mask_.MaskString());

  if (c4fileIn.empty())
    setupForWaterModel( wmIn );
  else {
    if (read_c4(c4fileIn)) return 1;
  }

  // Sanity check
  if (c4params_.empty()) {
    mprinterr("Error: No C4 parameters loaded.\n");
    return 1;
  }

  std::string polfile;
  if (polfileIn.empty()) {
    // First, need to determine where the Amber FF files are
    const char* env = getenv("AMBERHOME");
    if (env != 0)
      polfile = std::string(env) + "/dat/leap/parm/lj_1264_pol.dat";
    else {
      mprinterr("Error: AMBERHOME not set. Cannot find 'lj_1264_pol.dat'.\n");
      return 1;
    }
  } else
    polfile = polfileIn;
  if (!File::Exists(polfile)) {
    mprinterr("Error: Polarizability file not found: %s\n", polfile.c_str());
    return 1;
  }
  if (read_pol(polfile)) return 1;

  return 0;
}
/*
int LJ1264_Params::LoadNonbonds(std::string const& frcmodNameIn) const {
  std::string frcmodName;
  if (frcmodNameIn.empty()) {
    // Determine frcmod file from the water model
    // First, need to determine where the Amber FF files are
    const char* env = getenv("AMBERHOME");
    if (env != 0)
      frcmodName = std::string(env) + "/dat/leap/parm/frcmod.";
    else {
      mprinterr("Error: AMBERHOME not set. Cannot find LJ 12-6-4 frcmod files.\n");
      return 1;
    }

    switch (waterModel_) {
      case UNKNOWN_WATER_MODEL :
        mprintf("Warning: Unknown water model; cannot determine LJ 12-6-4 frcmod file named.\n");
        return 0;
      case Cpptraj::Parm::TIP3P : frcmodName.append("ions234lm_1264_tip3p");
*/

/** \return Help text. */
std::string LJ1264_Params::HelpText() {
  std::string out("[lj1264mask <mask>] [tunfactor <fac>] [c4file <c4file>] [polfile <polfile>] [lj1264solvent <solvent>]");
  return out;
}

/** Init from argument list */
int LJ1264_Params::Init_LJ1264(ArgList& argIn, int debugIn, std::string const& solventBoxNameIn) {
  std::string lj1264mask = argIn.GetStringKey("lj1264mask");
  std::string c4file = argIn.GetStringKey("c4file");
  std::string polfile = argIn.GetStringKey("polfile");
  double tunfactor = argIn.getKeyDouble("tunfactor", 1.0);
  std::string solventBoxName;
  if (solventBoxNameIn.empty())
    solventBoxName = argIn.GetStringKey("lj1264solvent");
  else
    solventBoxName = solventBoxNameIn;
  // Try to guess the water model if we are solvating
  Cpptraj::Parm::WaterModelType wm = Cpptraj::Parm::UNKNOWN_WATER_MODEL;
  if (!solventBoxName.empty()) {
    if (solventBoxName == "TIP3PBOX") wm = Cpptraj::Parm::TIP3P;
    else if (solventBoxName == "TIP4PEWBOX") wm = Cpptraj::Parm::TIP4PEW;
    else if (solventBoxName == "SPCBOX") wm = Cpptraj::Parm::SPCE;
    else if (solventBoxName == "OPC3BOX") wm = Cpptraj::Parm::OPC3;
    else if (solventBoxName == "OPCBOX") wm = Cpptraj::Parm::OPC;
    else if (solventBoxName == "FB3BOX") wm = Cpptraj::Parm::FB3;
    else if (solventBoxName == "FB4BOX") wm = Cpptraj::Parm::FB4;
    else
      mprintf("Warning: Unable to determine solvent model for LJ 12-6-4 from solvent box name %s\n",
              solventBoxName.c_str());
  }
  if (wm == Cpptraj::Parm::UNKNOWN_WATER_MODEL) {
    mprintf("Warning: Unable to determine water model for LJ 12-6-4. Using TIP3P.\n");
    wm = Cpptraj::Parm::TIP3P;
  }
  if (Init_LJ1264(lj1264mask, c4file, wm, polfile, tunfactor, debugIn)) {
    mprinterr("Error: Init of LJ 12-6-4 failed.\n");
    return 1;
  }
  return 0;
}

/** Parameterize */
int LJ1264_Params::AssignLJ1264(Topology& topOut)
{
  mprintf("\tAssigning LJ 12-6-4 parameters.\n");
  // Must already have nonbond info
  if (!topOut.Nonbond().HasNonbond()) {
    mprinterr("Error: Topology %s does not have LJ A/B nonbond parameters; required before setting LJ 12-6-4 nonbond parameters.\n");
    return 1;
  }
  if (topOut.SetupIntegerMask( mask_ )) {
    mprinterr("Error: Could not set up LJ 12-6-4 mask '%s'\n", mask_.MaskString());
    return 1;
  }
  mprintf("\t  LJ 12-6-4 mask selection: ");
  mask_.BriefMaskInfo();
  mprintf("\n");
  if (mask_.None()) {
    mprintf("Warning: Nothing selected by LJ 12-6-4 mask '%s'\n", mask_.MaskString());
    return 0;
  }

  typedef std::pair<int, double> IdxParamMapPair;
  typedef std::map<int, double> IdxParamMapType;
  //typedef std::pair<IdxParamMapType::iterator, bool> IdxParamMapRet;
  IdxParamMapType idxC4Map;

  // For each selected atom get the type index, element, charge, C4 parameter
  for (AtomMask::const_iterator at = mask_.begin(); at != mask_.end(); ++at)
  {
    int typeIdx = topOut[*at].TypeIndex();
    // Sanity check
    if (typeIdx < 0) {
      mprinterr("Error: Atom %s does not have an atom type index.\n", topOut.AtomMaskName(*at).c_str());
      return 1;
    }
    // Was this type index already seen?
    IdxParamMapType::const_iterator ret = idxC4Map.find( typeIdx );
    if (ret != idxC4Map.end()) {
      continue; // TODO check that C4 param matches?
    }
    // Figure out charge and element
    int iCharge = (int)topOut[*at].Charge();
    std::string eltName = std::string(topOut[*at].ElementName());
    // LJ 12-6-4 convention is element name second char is lower case
    if (eltName.size() > 1) eltName[1] = tolower(eltName[1]);
    // Append the integer charge
    std::string ename = eltName + integerToString(iCharge);
    // Get the C4 parameter
    NameMapType::const_iterator c4it = c4params_.find( ename );
    if (c4it == c4params_.end()) {
      mprinterr("Error: Could not find a C4 parameter for atom %s, type index %i, charge %i, element %s\n",
                topOut.AtomMaskName(*at).c_str(), typeIdx, iCharge, eltName.c_str());
      return 1;
    }
    // Save the type index/C4 param pair
    idxC4Map.insert( IdxParamMapPair(typeIdx, c4it->second) );

    if (debug_ > 0)
       mprintf("DEBUG: %s : idx=%i  chg=%i  ename=%s  C4=%g\n", topOut.AtomMaskName(*at).c_str(),
               typeIdx, iCharge, ename.c_str(), c4it->second);
  }
  mprintf("\t  %zu unique C4 parameters for atoms selected by '%s'\n", idxC4Map.size(), mask_.MaskString());

  // For each atom type index find the polarizability
  //IdxParamMapType idxPolMap;
  class TypePol {
    public:
      TypePol() : pol_(0) {}

      bool IsSet() const { return typeName_.len() > 0; }

      NameType typeName_;
      double pol_;
  };
  typedef std::vector<TypePol> IdxPolArray;
  IdxPolArray idxPol( topOut.Nonbond().Ntypes() );

  for (int at = 0; at != topOut.Natom(); at++)
  {
    Atom const& ATM = topOut[at];
    int typeIdx = ATM.TypeIndex();
    // Sanity check
    if (typeIdx < 0) {
      mprinterr("Error: Atom %s does not have an atom type index.\n", topOut.AtomMaskName(at).c_str());
      return 1;
    }
    // Find a polarization from the atom type
    NameMapType::const_iterator polit = pol_.find( ATM.Type().Truncated() );
    if (polit == pol_.end()) {
      mprinterr("Error: Could not find LJ 12-6-4 polarization for atom type %s (%s)\n", *(ATM.Type()), topOut.AtomMaskName(at).c_str());
      return 1;
    }
    // See if this type index was already seen
    if (idxPol[typeIdx].IsSet()) {

//    IdxParamMapType::const_iterator ret = idxPolMap.find( typeIdx );
//    if (ret != idxPolMap.end()) {
      // Ensure polarization is a match
      if ( fabs( polit->second - idxPol[typeIdx].pol_ ) > Constants::SMALL ) {
        mprinterr("Error: Polarizability parameter of atom type %s (%g) is not the same as that of type %s (%g),\n",
                  *(ATM.Type()), polit->second, *(idxPol[typeIdx].typeName_), idxPol[typeIdx].pol_);
        mprinterr("Error:   but their VDW parameters are the same.\n");
        return 1;
      }
    } else {
      // New polarization
      if (debug_ > 0)
        mprintf("DEBUG: idx= %i  type= %s  pol= %g\n", typeIdx, *(ATM.Type()), polit->second);
      idxPol[typeIdx].typeName_ = ATM.Type();
      idxPol[typeIdx].pol_ = polit->second;
      //idxPolMap.insert( IdxParamMapPair(typeIdx, polit->second) );
    }
  }

  // Allocate space for the LJC params
  topOut.SetNonbond().AllocateLJC();

  // Loop over all metal type indices
  for (IdxParamMapType::const_iterator it = idxC4Map.begin(); it != idxC4Map.end(); ++it)
  {
    // Loop over all other atom type indices
    for (int jidx = 0; jidx < topOut.Nonbond().Ntypes(); jidx++)
    {
      int ljidx = topOut.Nonbond().GetLJindex( it->first, jidx );
      //IdxParamMapType::const_iterator jt = idxPolMap.find( jidx );
      TypePol const& POL = idxPol[jidx];
      if (debug_ > 0)
        mprintf("DEBUG: Type1= %i  C4= %g  Type2= %i  NBidx= %i  Pol= %g\n", it->first, it->second,
                jidx, ljidx, POL.pol_);
      if (ljidx < 0) {
        mprinterr("Error: Type pair %i -- %i is an LJ 10-12 pair, not LJ 6-12.\n", it->first, jidx);
        return 1;
      }
      // Calculate final C4 term
      double ljc;
      if (POL.typeName_ == "OW")
        // There is only one C4 term that exists between water and a certain ion
        ljc = it->second;
      else
        // There are two C4 terms that need to add together between two different ions
        ljc = topOut.Nonbond().LJC_Array( ljidx ) + (it->second / WATER_POL_ * POL.pol_ * tunfactor_);
      if (debug_ > 0)
        mprintf("DEBUG:\t\tFinal LJC = %g\n", ljc);
      topOut.SetNonbond().SetLJC( ljidx, ljc );
    } // END loop over all atom types indices
  } // END loop over all metal type indices

  return 0;
}
