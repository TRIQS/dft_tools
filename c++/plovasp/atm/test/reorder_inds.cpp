
#include "testing.hpp"

int main()
{
  double e[4], en;
  int inds[4], inds_should[4];
  int flag;
  
  e[0] = -1.5;
  e[1] = -1.0;
  e[2] = -1.309017;
  e[3] = -0.5;

  en = -0.55;
  printf("\n  Test case 2\n\n");

  flag = dos_reorder(en, e, inds);

  inds_should[0] = 0;
  inds_should[1] = 2;
  inds_should[2] = 1;
  inds_should[3] = 3;

  if(check_reorder_inds(inds, inds_should)) return 1;
}

