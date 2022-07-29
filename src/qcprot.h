#ifndef INC_QCPROT_H
#define INC_QCPROT_H
namespace Cpptraj {
namespace QCPRot {
/*******************************************************************************
 *  -/_|:|_|_\-
 *
 *  File:           qcprot.h
 *  Version:        1.5
 *
 *  Function:       Rapid calculation of the least-squares rotation using a
 *                  quaternion-based characteristic polynomial and
 *                  a cofactor matrix
 *
 *  Author(s):      Douglas L. Theobald
 *                  Department of Biochemistry
 *                  MS 009
 *                  Brandeis University
 *                  415 South St
 *                  Waltham, MA  02453
 *                  USA
 *
 *                  dtheobald@brandeis.edu
 *
 *                  Pu Liu
 *                  Johnson & Johnson Pharmaceutical Research and Development, L.L.C.
 *                  665 Stockton Drive
 *                  Exton, PA  19341
 *                  USA
 *
 *                  pliu24@its.jnj.com
 *
 *
 *    If you use this QCP rotation calculation method in a publication, please
 *    reference:
 *
 *      Douglas L. Theobald (2005)
 *      "Rapid calculation of RMSD using a quaternion-based characteristic
 *      polynomial."
 *      Acta Crystallographica A 61(4):478-480.
 *
 *      Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2009)
 *      "Fast determination of the optimal rotational matrix for macromolecular
 *      superpositions."
 *      Journal of Computational Chemistry 31(7):1561-1563.
 *
 *  Copyright (c) 2009-2016 Pu Liu and Douglas L. Theobald
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without modification, are permitted
 *  provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this list of
 *    conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice, this list
 *    of conditions and the following disclaimer in the documentation and/or other materials
 *    provided with the distribution.
 *  * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to
 *    endorse or promote products derived from this software without specific prior written
 *    permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *  Source:         started anew.
 *
 *  Change History:
 *    2009/04/13      Started source
 *
 ******************************************************************************/

/* Calculate the RMSD & rotational matrix.

        Input:
               coords1 -- reference structure
               coords2 -- candidate structure
               len     -- the size of the system
               weight  -- the weight array of size len; set to NULL if not needed
        Output:
               rot[9]  -- rotation matrix
        Return:
               RMSD value

*/
//double CalcRMSDRotationalMatrix(double **coords1, double **coords2, const int len, double *rot, const double *weight);

/* Calculate the inner product of two structures.
   If weight array is not NULL, calculate the weighted inner product.

        Input:
               coords1 -- reference structure
               coords2 -- candidate structure
               len     -- the size of the system
               weight  -- the weight array of size len: set to NULL if not needed
        Output:
               A[9]    -- the inner product matrix
        Return:
                (G1 + G2) * 0.5; used as E0 in function 'FastCalcRMSDAndRotation'

        Warning:
            1. You MUST center the structures, coords1 and coords2, before calling this function.

            2. Please note how the structure coordinates are stored in the double **coords
               arrays. They are 3xN arrays, not Nx3 arrays as is also commonly
               used (where the x, y, z axes are interleaved). The difference is
               something like this for storage of a structure with 8 atoms:

           Nx3: xyzxyzxyzxyzxyzxyzxyzxyz
           3xN: xxxxxxxxyyyyyyyyzzzzzzzz

           The functions can be easily modified, however, to accomodate any
           data format preference. I chose this format because it is readily
           used in vectorized functions (SIMD, Altivec, MMX, SSE2, etc.).
*/
//double InnerProduct(double *A, double **coords1, double **coords2, const int len, const double *weight);

/** Calculate the RMSD, and/or the optimal rotation matrix.
  * \param rot  -- Output rotation matrix (order of xx, xy, xz, yx, yy, yz, zx, zy, zz).
  *                If null, do not calculate rotation matrix.
  * \param A[9] -- the inner product (coordinate covariance matrix) of two structures.
  * \param rmsd -- Output RMSD value.
  * \param E0   -- (G1 + G2) * 0.5
  * \param len  -- The size of the system (total mass if mass-weighted, # atoms otherwise).
  * \param minScore -- if( minScore > 0 && rmsd < minScore) then calculate only the rmsd;
                           otherwise, calculate both the RMSD & the rotation matrix.
  * \return -1 if only the RMSD was calculated.
  * \return >= 0 If rotation matrix and RMSD were calculated.
  */
int FastCalcRMSDAndRotation(double *rot, double *A, double *rmsd, double E0, double len);

/* Center the coordinates.

        Warning:
            If you are doing a full superposition (the usual least squares way),
            you MUST center each structure first. That is, you must translate
            each structure so that its centroid is at the origin.
            You can use CenterCoords() for this.
*/
//void CenterCoords(double **coords, const int len);
}
}
#endif
