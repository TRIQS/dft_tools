
#include "testing.hpp"

int main()
{
  double e[4], en, ci_sum, ct, res[4];
  int inds[4], inds_should[4];
  int i, flag;
  char mess[128];
  
  e[0] = -1.5;
  e[1] = -1.309017;
  e[2] = -1.0;
  e[3] = -0.5;

  en = -0.55;
  printf("\n  Test case 2\n\n");

  flag = dos_reorder(en, e, inds);

  dos_corner_weights(en, e, inds, res);
  dos_tet_weights(en, e, inds, &ct);

  for(i = 0, ci_sum = 0.0; i < 4; i++)
  {
    printf(" res[%d] = %20.15lf\n", i, res[i]);
    ci_sum += res[i];
  }

  printf("  Difference: %le\n", fabs(ci_sum - ct));

  sprintf(mess, "Difference between 'ci_sum' and 'ct' is too large");
  ASSERT(fabs(ci_sum - ct) < 1e-12, mess);
}

