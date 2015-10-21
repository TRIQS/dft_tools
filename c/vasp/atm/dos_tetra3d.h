#ifndef __C_DOS_TETRA3D_H__
#define __C_DOS_TETRA3D_H__

#include <Python.h>
#include <arrayobject.h>

static PyObject *tetra_DOS3D(PyObject *self, PyObject *args);

void tet_dos3d(double en, double *eigk, int strd_eigk,
                 npy_int64 *itt, int ntet, int *strd_itt,
                 double *cti, int *strd_cti);
int dos_corner_weights(double en, double *eigs, int *inds,
                 double *ci);

int dos_reorder(double en, double *e, int *inds);


#endif
