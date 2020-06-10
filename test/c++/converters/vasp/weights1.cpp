
#include "testing.hpp"

int main()
{
  double e[4], en, res[4], r_should[4];
  int inds[4];
  int flag;
  
  e[0] = -1.5;
  e[1] = -1.309017;
  e[2] = -1.0;
  e[3] = -0.5;

  en = -0.55;
  printf("\n  Test case 4\n\n");

  dos_corner_weights(en, e, inds, res);

  r_should[0] = 0.000309016992226;
  r_should[1] = 0.000381966005939;
  r_should[2] = 0.000618033984453;
  r_should[3] = 0.017232002550965;

  if(check_weights_result(res, r_should)) return 1;
}

