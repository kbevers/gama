/*
    C++ Matrix/Vector templates (gMatVec 0.9.14)
    Copyright (C) 2000  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the gMatVec C++ Matrix/Vector template library.

    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
 *  $Id: hilbert.h,v 1.1 2001/12/07 11:59:45 cepek Exp $
 */

#ifndef gMatVec__Hilbert_Matrix__h__
#define gMatVec__Hilbert_Matrix__h__

#include <gmatvec/gmatvec.h>

namespace gMatVec {

template <class Float, class Exc> Mat<Float, Exc> InvHilbert(Index n)
{
  /* Inverse of a finite segment of the Hilbert matrix
   *
   * CACM Algorithm 50 by John R. Herndon (1961), corrections by
   * B. Randel, improved version by P. Naur
   * ----------------------------------------------------------------- 
   *
   * This procedure computes the elements of the innverse of an n x n
   * finite segment of the Hilbert Matrix. The Hilbert matrix has the
   * elements H(i,j) = 1/(i+j-1). The segments of this are known to be
   * increasingly ill-conditioned with increasing size.
   * 
   */

  Mat<Float, Exc> H(n,n);

  Index i, j, k;
  Float w;

  w = H(1,1) = n*n;

  for (i=2; i<=n; i++)
    {
      Float k = Float(n+i-1)*(n-i+1)/((i-1)*(i-1));
      w = H(i,i) = w * k * k;
    }

  for (i=1; i<=n-1; i++)
    for (j=i+1; j<=n; j++)
      {
        k = j - 1;
        H(i,j) = -H(i,k)*(n+k)*(n-k)/(k*k);
      }

  for (i=2; i<=n; i++)
    for (j=1; j<=i; j++)
      H(i,j) = H(j,i) = H(j,i)/(i+j-1);

  return H;
}

template <class Float, class Exc> Mat<Float, Exc> Hilbert(Index n)
{
  gMatVec::Mat<Float, Exc> H(n,n);
  const Float one = 1;

  for (Index i=1; i<=n; i++)
    for (Index j=1; j<=n; j++)
      H(i,j) = one / (i+j-1);

  return H;
}

}

#endif


