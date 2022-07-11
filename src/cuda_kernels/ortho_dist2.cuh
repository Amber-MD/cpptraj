#ifndef INC_ORTHO_DIST2_CUH
#define INC_ORTHO_DIST2_CUH
/** \return Shortest imaged distance^2 between given coordinates in an orthorhombic box. */
template <typename T>
__device__ T ortho_dist2(T a1x, T a1y, T a1z,
                         T a2x, T a2y, T a2z,
                         const T* box)
{
  T x = a1x - a2x;
  T y = a1y - a2y;
  T z = a1z - a2z;
  // Get rid of sign info
  if (x<0) x=-x;
  if (y<0) y=-y;
  if (z<0) z=-z;
  // Get rid of multiples of box lengths 
  while (x > box[0]) x = x - box[0];
  while (y > box[1]) y = y - box[1];
  while (z > box[2]) z = z - box[2];
  // Find shortest distance in periodic reference
  T D2 = box[0] - x;
  if (D2 < x) x = D2;
  D2 = box[1] - y;
  if (D2 < y) y = D2;
  D2 = box[2] - z;
  if (D2 < z) z = D2;

  return (x*x + y*y + z*z);
}
#endif
