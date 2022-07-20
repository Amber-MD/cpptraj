#ifndef INC_ENERGY_ENE_ANGLE_H
#define INC_ENERGY_ENE_ANGLE_H
namespace Cpptraj {
namespace Energy {

template <typename T>
T Ene_Angle(T const* xyz1, T const* xyz2, T const* xyz3, T const& teq, T const& tk)
{
  static const T pt999 = 0.9990;
  T ij[3], kj[3];

  ij[0] = xyz1[0] - xyz2[0];
  ij[1] = xyz1[1] - xyz2[1];
  ij[2] = xyz1[2] - xyz2[2];
  kj[0] = xyz3[0] - xyz2[0];
  kj[1] = xyz3[1] - xyz2[1];
  kj[2] = xyz3[2] - xyz2[2];

  T rij = ij[0]*ij[0] + ij[1]*ij[1] + ij[2]*ij[2];
  T rkj = kj[0]*kj[0] + kj[1]*kj[1] + kj[2]*kj[2];

  T rik = sqrt(rij*rkj);
  T ct0 = (ij[0]*kj[0] + ij[1]*kj[1] + ij[2]*kj[2]) / rik;
  if (ct0 < -pt999)
    ct0 = -pt999;
  else if (ct0 > pt999)
    ct0 = pt999;
  //T ct1 = std::max(-pt999, ct0);
  //T ct2 = std::min( pt999, ct1);

  //T cst = ct2;
  T theta = acos(ct0);

  T da = theta - teq;
  // for rms deviation from ideal angles:
  // eadev = eadev + da*da
  T df = tk * da;
  T eaw = df * da;
  return eaw;
}
}
}
#endif
