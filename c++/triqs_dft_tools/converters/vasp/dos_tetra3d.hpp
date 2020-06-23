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
#pragma once

#include <triqs/arrays.hpp>


/// DOS of a band by analytical tetrahedron method
///
///   Returns corner weights for all tetrahedra for a given band and real energy.
triqs::arrays::array<double, 2>
dos_tetra_weights_3d(triqs::arrays::array_view<double, 1> eigk, /// Band energies for each k-point
                     double en, /// Energy at which DOS weights are to be calculated
                     triqs::arrays::array_view<long, 2> itt /// Tetrahedra defined by k-point indices
);
//array<double, 2> 
//dos_tetra_weights_3d(array<double, 1> eigk, /// Band energies for each k-point
//                     double e, /// Energy at which DOS weights are to be calculated
//                     array<long, 2> itt /// Tetrahedra defined by k-point indices
//);
int dos_corner_weights(double en, double *eigs, int *inds, double *ci);
int dos_tet_weights(double en, double *eigs, int *inds, double *ct);
int dos_reorder(double en, double *e, int *inds);



