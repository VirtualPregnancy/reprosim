
%module(package="reprosim") geometry
%include symbol_export.h

%{
#include "geometry.h"
%}

// define_rad_from_file has an optional argument that C cannot replicate,
// so we use SWIG to override with a C++ version that can.
void add_matching_mesh(const char *umbilical_elem_option, const char *UMB_ELEMS_FILE="");
void define_rad_from_file(const char *FIELDFILE, const char *order_system="strahler", double s_ratio=1.54);
void define_rad_from_geom(const char *ORDER_SYSTEM, double CONTROL_PARAM, const char *START_FROM, double START_RAD, const char *group_type_in="all", const char *group_option_in="dummy");

%include geometry.h

