#include <stdlib.h>

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


