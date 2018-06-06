#ifndef REPROSIM_PRESSURE_RESISTANCE_FLOW_H
#define REPROSIM_PRESSURE_RESISTANCE_FLOW_H

#include "symbol_export.h"

SHO_PUBLIC void evaluate_prq(const char *mesh_type, const char *bc_type, double inlet_flow, double inlet_pressure, double outlet_pressure);


#endif /* REPROSIM_PRESSURE_RESISTANCE_FLOW_H */
