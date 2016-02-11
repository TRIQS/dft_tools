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
#ifndef __C_DOS_TETRA3D_H__
#define __C_DOS_TETRA3D_H__

#include <Python.h>
#include <numpy/arrayobject.h>

static PyObject *tetra_DOS3D(PyObject *self, PyObject *args);

//void tet_dos3d(double en, double *eigk, int strd_eigk,
//                 npy_int64 *itt, int ntet, int *strd_itt,
//                 double *cti, int *strd_cti);
void tet_dos3d(double en, PyArrayObject *py_eigk,
                 PyArrayObject *py_itt, int ntet,
                 PyArrayObject *py_cti);
int dos_corner_weights(double en, double *eigs, int *inds,
                 double *ci);

int dos_reorder(double en, double *e, int *inds);

static const double small = 2.5e-2, tol = 1e-8;

#endif
