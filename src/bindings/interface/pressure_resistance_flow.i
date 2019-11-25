
%module(package="reprosim") pressure_resistance_flow
%include symbol_export.h
%include pressure_resistance_flow.h

%{
#include "pressure_resistance_flow.h"
%}

// some functions have an optional argument that C cannot replicate,
// so we use SWIG to override with a C++ version that can.
void calculate_stats(const char *FLOW_GEN_FILE, double image_voxel_size=0.00);
