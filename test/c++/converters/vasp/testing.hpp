
#include <cstdio>
#include <cstring>
#include <cmath>

#define ASSERT(cond, message) if(!(cond)) { \
        printf("*** Fail: %s\n", message); return 1;}

//
// Functions defined in 'atm_c' library
//
extern int dos_corner_weights(double en, double *eigs, int *inds,
                   double *ci);
extern int dos_tet_weights(double en, double *eigs, int *inds,
                    double *ct);
extern int dos_reorder(double en, double *e, int *inds);


int check_reorder_flag(int flag, int flag_should);
int check_reorder_inds(int *inds, int *inds_should);
int check_weights_result(double *res, double *r_should);

//
//  Test templates
//
int check_reorder_flag(int flag, int flag_should)
{
  char mess[128];
  sprintf(mess, "Reorder flag should be %d and not %d", flag_should, flag);
  ASSERT(flag == flag_should, mess);
  return 0;
}

int check_reorder_inds(int *inds, int *inds_should)
{
  char mess[128], tmp[128], numb[128];
  int i, flag;

  strcpy(mess, "Inds should be ");
  numb[0] = '\0';
  flag = 1;
  for(i = 0; i < 4; i++)
  {
    if(inds_should[i] != inds[i]) flag = 0;
    sprintf(tmp, "  %d", inds_should[i]);
    strcat(numb, tmp);
  }
  strcat(mess, numb);
  strcat(mess, " and not");
  numb[0] = '\0';
  for(i = 0; i < 4; i++)
  {
    sprintf(tmp, "  %d", inds[i]);
    strcat(numb, tmp);
  }
  strcat(mess, numb);
  
  ASSERT(flag, mess);

  return 0;
}

int check_weights_result(double *res, double *r_should)
{
  const double tol = 1e-14;
  int i, flag;
  char mess[128], tmp[128];

  flag = 1;
  for(i = 0; i < 4; i++)
  {
    if(fabs(r_should[i] - res[i]) > tol)
    {
      flag = 0;
      break;
    }
  }
  strcpy(mess, "Success");
  if(!flag)
  {
    sprintf(mess, "res[%d] should be %20.15lf", i, r_should[i]);
    sprintf(tmp, " and not %20.15lf", res[i]);
    strcat(mess, tmp);
  }
  
  ASSERT(flag, mess);

  return 0;
}


