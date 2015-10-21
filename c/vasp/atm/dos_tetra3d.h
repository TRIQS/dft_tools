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
