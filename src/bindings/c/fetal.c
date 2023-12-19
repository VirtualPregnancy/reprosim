#include "fetal.h"
#include "utils.h"
#include <string.h>


void fetal_model_c(const char *OUTDIR, int *filename_len, double *dt,int *num_heart_beats,double *T_beat,double *T_vs,double *T_as,double *T_v_delay,double *U0RV,double *EsysRV,double *EdiaRV,double *RvRv,double *U0LV,double *EsysLV,double *EdiaLV,double *RvLV,double *U0A, double *V0V, double *V0A);
void assign_fetal_arrays_c();

void fetal_model(const char *OUTDIR, double dt,int num_heart_beats,double T_beat,double T_vs,double T_as,double T_v_delay,double U0RV,double EsysRV,double EdiaRV,double RvRv,double U0LV,double EsysLV,double EdiaLV,double RvLV,double U0A, double V0V, double V0A)
{
  int filename_len = strlen(OUTDIR);
  fetal_model_c(OUTDIR, &filename_len,&dt,&num_heart_beats,&T_beat,&T_vs,&T_as,&T_v_delay,&U0RV,&EsysRV,&EdiaRV,&RvRv,&U0LV,&EsysLV,&EdiaLV,&RvLV,&U0A,&V0V,&V0A);
}

void assign_fetal_arrays()
{
  assign_fetal_arrays_c();
}