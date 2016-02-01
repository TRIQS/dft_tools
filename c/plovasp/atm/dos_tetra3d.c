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
#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <numpy/arrayobject.h>
#include <complex.h>
//#include "tetra.h"
#include "argsort.h"
#include "dos_tetra3d.h"

/***************************************************

  Analytical tetrahedron method as described in
Lambin et al., PRB 29, 6, 3430 (1984).

***************************************************/

static double F(double en, double e1, double e2, double e3, double e4);
static double K2(double en, double e1, double e2, double e3);
static double K1(double en, double e1, double e2);

static void fun_dos_case1(double en, double *eigs, double *ci);
static void fun_dos_case2(double en, double *eigs, double *ci);
static void fun_dos_case3(double en, double *eigs, double *ci);

int dos_corner_weights(double en, double *eigs, int *inds, double *ci);
int dos_tet_weights(double en, double *eigs, int *inds, double *ct);

const static int NUM_TET_CORNERS = 4;

static PyMethodDef c_tetraMethods[] = {
  {"dos_weights_3d", tetra_DOS3D, METH_VARARGS,
   "C-implementation of the tetrahedron method for calculating DOS"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initc_atm_dos(void)
{
  (void) Py_InitModule("c_atm_dos", c_tetraMethods);
  import_array();
}

//
// Integration_Weights
//
static PyObject *
tetra_DOS3D(PyObject *self, PyObject *args)
{

  PyArrayObject *py_eigk, *py_itt;  // Input Numpy arrays
  PyArrayObject *py_cti;  // Output Numpy array
//  npy_int64 *itt;  // C-pointer to the 'itt' array
  npy_intp *dims;  // Dimensions of the output array 'rti'
  double *cti; // C-pointer to 'eigk' and 'rti' arrays
  double e;
//  double et[4]; // Real part of eigenvalues
//  int strd_itt[2], strd_eigk, strd_cti[2], it_cti; // Strides
  int ntet; // Indices of sorted 'eigs'; number of tetrahedra
// Auxiliary variables and loop indices
  int nd;

// Data-preparation part
  if (!PyArg_ParseTuple(args, "O!dO!", 
      &PyArray_Type, &py_eigk, &e, &PyArray_Type, &py_itt))
    return NULL;

// Sanity tests (the types are assumed to be checked in the python wrapper)
//  nd = py_eigk->nd;
  nd = PyArray_NDIM(py_eigk);
  if (nd != 1)
  {
    PyErr_SetString(PyExc_ValueError, "  Array 'eigk' must be 1D");
    return NULL;
  }

//  nd = py_itt->nd;
  nd = PyArray_NDIM(py_itt);
  if (nd != 2)
  {
    PyErr_SetString(PyExc_ValueError, "  Array 'itt' must be 2D");
    return NULL;
  }

//  nd = py_itt->dimensions[0];
  dims = PyArray_DIMS(py_itt);
  if (dims[0] != 5)
  {
    PyErr_SetString(PyExc_ValueError, 
       "  The first dimension of 'itt' must be equal to 5");
    return NULL;
  }

//  ntet = (int) py_itt->dimensions[1];
  ntet = dims[1];

//  eigk = (double *)py_eigk->data;
//  strd_eigk = py_eigk->strides[0] / sizeof(npy_float64);

//  itt = (npy_int64 *)py_itt->data;
//  strd_itt[0] = py_itt->strides[0] / sizeof(npy_int64);
//  strd_itt[1] = py_itt->strides[1] / sizeof(npy_int64);

// Resulting array (the question is whether dims is copied or not?)
  cti = (double *)malloc(NUM_TET_CORNERS * ntet * sizeof(double));
  dims = (npy_intp *)malloc(2 * sizeof(npy_intp));

  dims[0] = NUM_TET_CORNERS;
  dims[1] = ntet;

  py_cti = (PyArrayObject *)PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, cti);

//  strd_cti[0] = py_cti->strides[0] / sizeof(double);
//  strd_cti[1] = py_cti->strides[1] / sizeof(double);

// Main part: now we can fill the 'cti' array and it will be returned
// by 'py_cti' as a numpy array

//
// Main function
//
//  tet_dos3d(e, eigk, strd_eigk, itt, ntet, strd_itt, cti, strd_cti);
  tet_dos3d(e, py_eigk, py_itt, ntet, py_cti);
 
  return PyArray_Return(py_cti);
}

//void tet_dos3d(double en, double *eigk, int strd_eigk,
//                 npy_int64 *itt, int ntet, int *strd_itt,
//                 double *cti, int *strd_cti)
void tet_dos3d(double en, PyArrayObject *py_eigk, 
                 PyArrayObject *py_itt, int ntet,
                 PyArrayObject *py_cti)
{
  double eigs[4], ci[4];
  int i, it, ik, inds[4], flag;
// **** DEBUG
  double ct, ci_sum;

// Loop over tetrahedra (triangles)
  for (it = 0; it < ntet; it++)
  {
//    it_cti = it * strd_cti[1];
// Sort eigenvalues and obtain indices of the sorted array
//   eigs: sorted eigenvalues
//   inds: index map
    for (i = 1; i < 5; i++) 
    {
//      ik = itt[i * strd_itt[0] + it * strd_itt[1]];
      ik = ((int *)PyArray_GETPTR2(py_itt, i, it))[0];
//      eigs[i - 1] = eigk[ik * strd_eigk];
      eigs[i - 1] = ((double *)PyArray_GETPTR1(py_eigk, ik))[0];
    }

// corner weights for a single tetrahedron
    dos_corner_weights(en, eigs, inds, ci);
    for(i = 0, ci_sum = 0.0; i < 4; i++)
      ci_sum += ci[i];

    flag = dos_tet_weights(en, eigs, inds, &ct);
    if(fabs(ct - ci_sum) > tol)
      {
        printf("  *** Error in weights: it = %d, flag = %d, en = %lf", it, flag, en);
        for(i = 0; i < 4; i++)
          printf(", e[%d] = %lf", i, eigs[i]);
        printf(",    c_diff = %le\n", fabs(ct - ci_sum));
        return;
      }  
//    printf("  it = %d, flag = %d", it, j);
//    for(i = 0; i < 4; i++)
//      printf(", e[%d] = %lf", i, eigs[i]);
//    printf(",    c_diff = %le\n", fabs(ct - ci_sum));

//    if(j < 4)
//    {
//      printf("  flag = %d, e = %lf", j, en);
//      for(i = 0; i < 4; i++)
//        printf(", eigs[%d] = %lf", i, eigs[i]);
//      printf("\n");
//      printf("  ci = ");
//      for(i = 0; i < 4; i++)
//        printf(", %lf", ci[i]);
//      printf("\n");
//    }

    for(i = 0; i < 4; i++)
    {
//      j = inds[i] * strd_cti[0] + it_cti;
//      cti[j] = ci[i];
      ((double *)PyArray_GETPTR2(py_cti, inds[i], it))[0] = ci[i];
    }

  }  // it = 1, ntet
}
 
int dos_corner_weights(double en, double *eigs, int *inds,
                   double *ci)
{
  int flag, i;
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

int dos_tet_weights(double en, double *eigs, int *inds,
                    double *ct)
{
  double e1, e2, e3, e4;
  double complex s;
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
    if(fabs(e2 - e1) > tol && fabs(e3 - e1) > tol && fabs(e4 - e1) > tol)
      *ct = 3.0 * (en - e1) * (en - e1) / ((e2 - e1) * (e3 - e1) * (e4 - e1));
    else
    {
      s = fmin(fabs(e1 - e2), fabs(e3 - e1));
      s = fmin(fabs(s), fabs(e4 - e1));
      s /= 100.0;
      s = fmax(s, 1.0e-20) * I;

      *ct = 3.0 * creal((en - e1 + s) * (en - e1 + s) / ((e2 - e1 + s) * (e3 - e1 + s) * (e4 - e1 + s)));
    }
 
    break;

// E2 <= E <= E3
  case 2:
    if(fabs(e4 - e2) > tol && fabs(e3 - e2) > tol && fabs(e4 - e1) > tol && fabs(e3 - e1) > tol)
      *ct = 3.0 * (
          (e3 - en) * (en - e2) / ((e4 - e2) * (e3 - e2) * (e3 - e1)) +
          (e4 - en) * (en - e1) / ((e4 - e1) * (e4 - e2) * (e3 - e1)));
    else
    {
      s = fmin(fabs(e3 - e2), fabs(e3 - e1));
      s = fmin(fabs(s), fabs(e4 - e1));
      s = fmin(fabs(s), fabs(e4 - e2));
      s /= 100.0;
      s = fmax(s, 1.0e-20) * I;

      *ct = 3.0 * creal((
          (e3 - en + s) * (en - e2 + s) / ((e4 - e2 + s) * (e3 - e2 + s) * (e3 - e1 + s)) +
          (e4 - en + s) * (en - e1 + s) / ((e4 - e1 + s) * (e4 - e2 + s) * (e3 - e1 + s))));
    }
    break;

// E3 <= E <= E4
  case 3:
    if(fabs(e4 - e2) > tol && fabs(e4 - e3) > tol && fabs(e4 - e1) > tol)
      *ct = 3.0 * (e4 - en) * (e4 - en) / ((e4 - e1) * (e4 - e2) * (e4 - e3));
    else
    {
      s = fmin(fabs(e4 - e2), fabs(e4 - e1));
      s = fmin(fabs(s), fabs(e4 - e3));
      s /= 100.0;
      s = fmax(s, 1.0e-20) * I;

      *ct = 3.0 * creal((e4 - en + s) * (e4 - en + s) / ((e4 - e1 + s) * (e4 - e2 + s) * (e4 - e3 + s)));
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

int dos_reorder(double en, double *e, int *inds)
{
  double *ptrs[4], e_tmp[4];
  int i;

  for(i = 0; i < 4; i++)
    e_tmp[i] = e[i];

  argsort(e_tmp, inds, ptrs, 4);
  
  for(i = 0; i < 4; i++)
    e[i] = e_tmp[inds[i]];
 
  if((e[0] <= en && en <= e[3]) && fabs(e[3] - e[0]) < tol) return 6;
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
//  int i;

//  if(fabs(eigs[1] - eigs[0]) < tol)
//  {
//    for(i = 0; i < 4; i++) ci[i] = 0.0;
//    return;
//  }

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

//  if(fabs(eigs[2] - eigs[1]) < tol)
//  {
//    for(i = 0; i < 4; i++) ci[i] = 0.0;
//    return;
//  }

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

//  if(fabs(eigs[3] - eigs[2]) < tol)
//  {
//    for(i = 0; i < 4; i++) ci[i] = 0.0;
//    return;
//  }

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
  double complex s;

  if(fabs(e1 - e3) > tol && fabs(e4 - e2) > tol)
    return (e1 - en) * (en - e2) / ((e1 - e3) * (e4 - e2));
  else
  {
    s = fmin(fabs(e3 - e1), fabs(e4 - e2));
    s /= 100.0;
    s = fmax(s, 1.0e-20) * I;
   
    return creal((e1 - en + s) * (en - e2 + s) / ((e1 - e3 + s) * (e4 - e2 + s)));
  }
}

static double K2(double en, double e1, double e2, double e3)
{
  double complex s;

  if(fabs(e1 - e3) > tol && fabs(e1 - e2) > tol)
    return (en - e1) / ((e2 - e1) * (e3 - e1));
  else
  {
    s = fmin(fabs(e3 - e1), fabs(e1 - e2));
    s /= 100.0;
    s = fmax(s, 1.0e-20) * I;
   
    return creal((en - e1 + s) / ((e2 - e1 + s) * (e3 - e1 + s)));
  }
}

static double K1(double en, double e1, double e2)
{
  double complex s;

  if(fabs(e1 - e2) > tol)
    return (e1 - en) / ((e2 - e1) * (e2 - e1));
  else
  {
    s = fabs(e1 - e2);
    s /= 100.0;
    s = fmax(s, 1.0e-20) * I;
   
    return creal((e1 - en + s) / ((e2 - e1 + s) * (e2 - e1 + s)));
  }
}

