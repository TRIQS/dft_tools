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
#include <cstdlib>

int cmp(const void *a, const void *b)
{
  return (**(const double **)a) < (**(const double **)b) ? -1 : 1;
}

int icmp(const void *a, const void *b)
{
  return (**(const int **)a) - (**(const int **)b);
}

void argsort(double *arr, int *inds, double **ptrs, const int n)
{
  int i;
  
  for (i=0; i < n; i++)
  {
   ptrs[i] = arr + i;
  }

  qsort(ptrs, n, sizeof(double *), cmp);

  for (i=0; i < n; i++)
  {
   inds[i] = (int)(ptrs[i] - arr);
  }
}

void iargsort(int *iarr, int *inds, int **ptrs, const int n)
{
  int i;
  
  for (i=0; i < n; i++)
  {
   ptrs[i] = iarr + i;
  }

  qsort(ptrs, n, sizeof(int *), icmp);

  for (i=0; i < n; i++)
  {
   inds[i] = (int)(ptrs[i] - iarr);
  }
}


