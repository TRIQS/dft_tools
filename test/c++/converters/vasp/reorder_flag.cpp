
#include "testing.hpp"

int main()
{
  double e[4], en;
  int inds[4];
  int flag, flag_should;
  
  e[0] = -1.5;
  e[1] = -1.309017;
  e[2] = -1.0;
  e[3] = -0.5;

  en = -0.55;
  printf("\n  Test case 1\n\n");

  flag = dos_reorder(en, e, inds);
  flag_should = 3;
 
  if(check_reorder_flag(flag, flag_should)) return 1;
}

