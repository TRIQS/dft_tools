/*******************************************************************************
*
* This file is part of the ATM library.
*
* Copyright (C) 2010 by O. E. Peil
*
* TRIQS is free software: you can redistribute it and/or modify it under the
* terms of the GNU General Public License as published by the Free Software
* Foundation, either version 3 of the License, or (at your option) any later
* version.
*
* TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public License along with
* TRIQS. If not, see <http://www.gnu.org/licenses/>.
*
*******************************************************************************/
#include <triqs/arrays.hpp>
#include <iostream>
#include <complex>

#include "argsort.hpp"
#include "dos_tetra3d.hpp"

//#define __TETRA_DEBUG
#define __TETRA_ARRAY_VIEW

using triqs::arrays::array;
using triqs::arrays::array_view;

/***************************************************

  Analytical tetrahedron method as described in
Lambin et al., PRB 29, 6, 3430 (1984).

***************************************************/

/// Main function
//#ifdef __TETRA_ARRAY_VIEW
//void tet_dos3d(double en, array_view<double, 1>& eigk,
//                 array_view<long, 2>& itt, int ntet,
//                 array<double, 2>& cti);
//#else
//void tet_dos3d(double en, array<double, 1>& eigk,
//                 array<long, 2>& itt, int ntet,
//                 array<double, 2>& cti);
//#endif


/// Internal functions
int dos_corner_weights(double en, double *eigs, int *inds, double *ci);
int dos_tet_weights(double en, double *eigs, int *inds, double *ct);
int dos_reorder(double en, double *e, int *inds);

static double F(double en, double e1, double e2, double e3, double e4);
static double K2(double en, double e1, double e2, double e3);
static double K1(double en, double e1, double e2);

static void fun_dos_case1(double en, double *eigs, double *ci);
static void fun_dos_case2(double en, double *eigs, double *ci);
static void fun_dos_case3(double en, double *eigs, double *ci);

static const int NUM_TET_CORNERS = 4;
static const std::complex<double> I(0.0, 1.0);
static const double small = 2.5e-2, tol = 1e-8;

/*
  Returns corner contributions to the DOS of a band
*/
#ifdef __TETRA_ARRAY_VIEW
array<double, 2> dos_tetra_weights_3d(array_view<double, 1> eigk, double en, array_view<long, 2> itt)
#else
array<double, 2> dos_tetra_weights_3d(array<double, 1> eigk, double en, array<long, 2> itt)
#endif
{
  int ntet; /// Number of tetrahedra
// Auxiliary variables and loop indices

  if (first_dim(itt) != NUM_TET_CORNERS + 1)
  {
      TRIQS_RUNTIME_ERROR << "  The first dimension of 'itt' must be equal to 5";
  }

  ntet = second_dim(itt);

  array<double, 2> cti(NUM_TET_CORNERS, ntet); // Corner weights to be returned

//  tet_dos3d(e, eigk, itt, ntet, cti);

//
// Main algorithm (transferred from 'tet_dos3d()')
//
  double eigs[4], ci[4];
  
  int i, it, ik, inds[4], flag;
#ifdef __TETRA_DEBUG
  double ct, ci_sum;
#endif

// Loop over tetrahedra (triangles)
  for (it = 0; it < ntet; it++)
  {
    for (i = 1; i < 5; i++) 
    {
      ik = itt(i, it);
      eigs[i - 1] = eigk(ik);
    }

// Corner weights for a single tetrahedron
    dos_corner_weights(en, eigs, inds, ci);

#ifdef __TETRA_DEBUG
    for(i = 0, ci_sum = 0.0; i < 4; i++)
      ci_sum += ci[i];

    flag = dos_tet_weights(en, eigs, inds, &ct);
    if(std::abs(ct - ci_sum) > tol)
      {
        std::cout << "  *** Error in weights: it = " << it <<" flag = " << flag << ", en = " << en;
        for(i = 0; i < 4; i++)
          std::cout << ", e[" << i << "] = " << eigs[i];
        std::cout << ",    c_diff = " << std::abs(ct - ci_sum) << std::endl;

        TRIQS_RUNTIME_ERROR << "  Failed consistency check";
      }  
#endif

  for(i = 0; i < 4; i++)
    {
      cti(inds[i], it) = ci[i];
    }

  }  // it = 1, ntet

  return array_view<double,2>(cti);
}

//#ifdef __TETRA_ARRAY_VIEW
//void tet_dos3d(double en, array_view<double, 1>& eigk,
//                 array_view<long, 2>& itt, int ntet,
//                 array<double, 2>& cti)
//#else
//void tet_dos3d(double en, array<double, 1>& eigk,
//                 array<long, 2>& itt, int ntet,
//                 array<double, 2>& cti)
//#endif
//{
//  double eigs[4], ci[4];
//  
//  int i, it, ik, inds[4], flag;
//#ifdef __TETRA_DEBUG
//  double ct, ci_sum;
//#endif
//
//// Loop over tetrahedra (triangles)
//  for (it = 0; it < ntet; it++)
//  {
//    for (i = 1; i < 5; i++) 
//    {
//      ik = itt(i, it);
//      eigs[i - 1] = eigk(ik);
//    }
//
//// Corner weights for a single tetrahedron
//    dos_corner_weights(en, eigs, inds, ci);
//
//#ifdef __TETRA_DEBUG
//    for(i = 0, ci_sum = 0.0; i < 4; i++)
//      ci_sum += ci[i];
//
//    flag = dos_tet_weights(en, eigs, inds, &ct);
//    if(std::abs(ct - ci_sum) > tol)
//      {
//        std::cout << "  *** Error in weights: it = " << it <<" flag = " << flag << ", en = " << en;
//        for(i = 0; i < 4; i++)
//          std::cout << ", e[" << i << "] = " << eigs[i];
//        std::cout << ",    c_diff = " << std::abs(ct - ci_sum) << std::endl;
//        return;
//      }  
//#endif
//
//  for(i = 0; i < 4; i++)
//    {
//      cti(inds[i], it) = ci[i];
//    }
//
//  }  // it = 1, ntet
//}
 
/// Corner contributions to DOS
int dos_corner_weights(double en, double *eigs, int *inds,
                   double *ci)
{
  int flag, i;
// Sort eigenvalues and obtain indices of the sorted array
//   eigs: sorted eigenvalues
//   inds: index map
  flag = dos_reorder(en, eigs, inds);

  switch(flag)
  {
// E1 <= E <= E2
  case 1:
    fun_dos_case1(en, eigs, ci);
    break;

// E2 <= E <= E3
  case 2:
    fun_dos_case2(en, eigs, ci);
    break;

// E3 <= E <= E4
  case 3:
    fun_dos_case3(en, eigs, ci);
    break;

// E < E1 || E4 < E
  case 4:
  case 5:
    for(i = 0; i < 4; i++) ci[i] = 0.0;
    break;

// E1 == E4 == E
  case 6:
    for(i = 0; i < 4; i++) ci[i] = 0.25;
    break;
  }

  return flag;
}

/// Total (tetrahedron) contribution to DOS.
/// Here, it is calculated directly using an analytical formula.
/// This is mainly needed for debugging.
int dos_tet_weights(double en, double *eigs, int *inds,
                    double *ct)
{
  double e1, e2, e3, e4;
  std::complex<double> s;
  int flag;
  flag = dos_reorder(en, eigs, inds);

  e1 = eigs[0];
  e2 = eigs[1];
  e3 = eigs[2];
  e4 = eigs[3];

  switch(flag)
  {
// E1 <= E <= E2
  case 1:
    if(std::abs(e2 - e1) > tol && std::abs(e3 - e1) > tol && std::abs(e4 - e1) > tol)
      *ct = 3.0 * (en - e1) * (en - e1) / ((e2 - e1) * (e3 - e1) * (e4 - e1));
    else
    {
      s = fmin(std::abs(e1 - e2), std::abs(e3 - e1));
      s = fmin(std::abs(s), std::abs(e4 - e1));
      s /= 100.0;
      s = fmax(std::abs(s), 1.0e-20) * I;

      *ct = 3.0 * std::real((en - e1 + s) * (en - e1 + s) / ((e2 - e1 + s) * (e3 - e1 + s) * (e4 - e1 + s)));
    }
 
    break;

// E2 <= E <= E3
  case 2:
    if(std::abs(e4 - e2) > tol && std::abs(e3 - e2) > tol && std::abs(e4 - e1) > tol && std::abs(e3 - e1) > tol)
      *ct = 3.0 * (
          (e3 - en) * (en - e2) / ((e4 - e2) * (e3 - e2) * (e3 - e1)) +
          (e4 - en) * (en - e1) / ((e4 - e1) * (e4 - e2) * (e3 - e1)));
    else
    {
      s = fmin(std::abs(e3 - e2), std::abs(e3 - e1));
      s = fmin(std::abs(s), std::abs(e4 - e1));
      s = fmin(std::abs(s), std::abs(e4 - e2));
      s /= 100.0;
      s = fmax(std::abs(s), 1.0e-20) * I;

      *ct = 3.0 * std::real((
          (e3 - en + s) * (en - e2 + s) / ((e4 - e2 + s) * (e3 - e2 + s) * (e3 - e1 + s)) +
          (e4 - en + s) * (en - e1 + s) / ((e4 - e1 + s) * (e4 - e2 + s) * (e3 - e1 + s))));
    }
    break;

// E3 <= E <= E4
  case 3:
    if(std::abs(e4 - e2) > tol && std::abs(e4 - e3) > tol && std::abs(e4 - e1) > tol)
      *ct = 3.0 * (e4 - en) * (e4 - en) / ((e4 - e1) * (e4 - e2) * (e4 - e3));
    else
    {
      s = fmin(std::abs(e4 - e2), std::abs(e4 - e1));
      s = fmin(std::abs(s), std::abs(e4 - e3));
      s /= 100.0;
      s = fmax(std::abs(s), 1.0e-20) * I;

      *ct = 3.0 * std::real((e4 - en + s) * (e4 - en + s) / ((e4 - e1 + s) * (e4 - e2 + s) * (e4 - e3 + s)));
    }

    break;

// E < E1 || E4 < E
  case 4:
  case 5:
    *ct = 0.0;
    break;

// E1 == E4 == E
  case 6:
    *ct = 1.0;
    break;
  }

  return flag;
}

/// Sorts eigenvalues and also determines eigenvalue degeneracies.
/// Returns a case number corresponding to a combination of degeneracies.
int dos_reorder(double en, double *e, int *inds)
{
  double *ptrs[4], e_tmp[4];
  int i;

  for(i = 0; i < 4; i++)
    e_tmp[i] = e[i];

  argsort(e_tmp, inds, ptrs, 4);
  
  for(i = 0; i < 4; i++)
    e[i] = e_tmp[inds[i]];
 
  if((e[0] <= en && en <= e[3]) && std::abs(e[3] - e[0]) < tol) return 6;
  if(e[0] <= en && en <= e[1]) return 1;
  if(e[1] <= en && en <= e[2]) return 2;
  if(e[2] <= en && en <= e[3]) return 3;
  if(en < e[0]) return 4;
  if(e[3] < en) return 5;
  return -1;
}

static void fun_dos_case1(double en, double *eigs, double *ci)
{
  double e1, e2, e3, e4;

  e1 = eigs[0];
  e2 = eigs[1];
  e3 = eigs[2];
  e4 = eigs[3];

  ci[0] = K2(en, e1, e2, e4) * F(en, e2, e1, e1, e3) +
          K2(en, e1, e2, e3) * F(en, e3, e1, e1, e4) +
          K2(en, e1, e3, e4) * F(en, e4, e1, e1, e2);

  ci[1] = -K1(en, e1, e2) * F(en, e1, e1, e3, e4);

  ci[2] = -K1(en, e1, e3) * F(en, e1, e1, e2, e4);

  ci[3] = -K1(en, e1, e4) * F(en, e1, e1, e2, e3);
}

static void fun_dos_case2(double en, double *eigs, double *ci)
{
  double e1, e2, e3, e4;

  e1 = eigs[0];
  e2 = eigs[1];
  e3 = eigs[2];
  e4 = eigs[3];

  ci[0] = 0.5 * (K1(en, e3, e1) * (
          F(en, e3, e2, e2, e4) +
          F(en, e4, e1, e2, e4) +
          F(en, e3, e1, e2, e4)) +
                 K1(en, e4, e1) * (
          F(en, e4, e1, e2, e3) +
          F(en, e4, e2, e2, e3) +
          F(en, e3, e1, e2, e3)));

  ci[1] = 0.5 * (K1(en, e3, e2) * (
          F(en, e3, e2, e1, e4) +
          F(en, e4, e2, e1, e4) +
          F(en, e3, e1, e1, e4)) +
                 K1(en, e4, e2) * (
          F(en, e3, e2, e1, e3) +
          F(en, e4, e1, e1, e3) +
          F(en, e4, e2, e1, e3)));

  ci[2] = 0.5 * (-K1(en, e2, e3) * (
          F(en, e3, e2, e1, e4) +
          F(en, e4, e2, e1, e4) +
          F(en, e3, e1, e1, e4)) -
                  K1(en, e1, e3) * (
          F(en, e3, e2, e2, e4) +
          F(en, e4, e1, e2, e4) +
          F(en, e3, e1, e2, e4)));

  ci[3] = 0.5 * (-K1(en, e2, e4) * (
          F(en, e3, e2, e1, e3) +
          F(en, e4, e1, e1, e3) +
          F(en, e4, e2, e1, e3)) -
                  K1(en, e1, e4) * (
          F(en, e4, e1, e2, e3) +
          F(en, e4, e2, e2, e3) +
          F(en, e3, e1, e2, e3)));
}

static void fun_dos_case3(double en, double *eigs, double *ci)
{
  double e1, e2, e3, e4;

  e1 = eigs[0];
  e2 = eigs[1];
  e3 = eigs[2];
  e4 = eigs[3];

  ci[0] = K1(en, e4, e1) * F(en, e4, e4, e2, e3);
  
  ci[1] = K1(en, e4, e2) * F(en, e4, e4, e1, e3);

  ci[2] = K1(en, e4, e3) * F(en, e4, e4, e1, e2);

  ci[3] = -K2(en, e4, e3, e1) * F(en, e4, e3, e2, e4) -
           K2(en, e4, e2, e3) * F(en, e4, e2, e1, e4) -
           K2(en, e4, e1, e2) * F(en, e4, e1, e3, e4);

}

static double F(double en, double e1, double e2, double e3, double e4)
{
  std::complex<double> s;

  if(std::abs(e1 - e3) > tol && std::abs(e4 - e2) > tol)
    return (e1 - en) * (en - e2) / ((e1 - e3) * (e4 - e2));
  else
  {
// Regularization to avoid division by zero
    s = fmin(std::abs(e3 - e1), std::abs(e4 - e2));
    s /= 100.0;
    s = fmax(std::abs(s), 1.0e-20) * I;
   
    return std::real((e1 - en + s) * (en - e2 + s) / ((e1 - e3 + s) * (e4 - e2 + s)));
  }
}

static double K2(double en, double e1, double e2, double e3)
{
  std::complex<double> s;

  if(std::abs(e1 - e3) > tol && std::abs(e1 - e2) > tol)
    return (en - e1) / ((e2 - e1) * (e3 - e1));
  else
  {
// Regularization to avoid division by zero
    s = fmin(std::abs(e3 - e1), std::abs(e1 - e2));
    s /= 100.0;
    s = fmax(std::abs(s), 1.0e-20) * I;
   
    return std::real((en - e1 + s) / ((e2 - e1 + s) * (e3 - e1 + s)));
  }
}

static double K1(double en, double e1, double e2)
{
  std::complex<double> s;

  if(std::abs(e1 - e2) > tol)
    return (e1 - en) / ((e2 - e1) * (e2 - e1));
  else
  {
// Regularization to avoid division by zero
    s = std::abs(e1 - e2);
    s /= 100.0;
    s = fmax(std::abs(s), 1.0e-20) * I;
   
    return std::real((e1 - en + s) / ((e2 - e1 + s) * (e2 - e1 + s)));
  }
}


