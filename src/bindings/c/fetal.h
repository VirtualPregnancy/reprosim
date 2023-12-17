#ifndef REPROSIM_FETAL_H
#define REPROSIM_FETAL_H

#include "symbol_export.h"
SHO_PUBLIC void fetal_model(const char *OUTDIR, double dt, int num_heart_beats,double T_beat,double T_vs,double T_as,double T_v_delay,double U0RV,double EsysRV,double EdiaRV,double RvRv,double U0LV,double EsysLV,double EdiaLV,double RvLV,double U0A,double V0V, double V0A);
SHO_PUBLIC void assign_fetal_arrays();

#endif /* REPROSIM_FETAL_H */