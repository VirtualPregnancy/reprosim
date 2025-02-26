
#ifndef REPROSIM_GEOMETRY_H
#define REPROSIM_GEOMETRY_H

#include "symbol_export.h"

SHO_PUBLIC void add_matching_mesh(const char *umbilical_elem_option, int umbilical_element_numbers_len, int umbilical_element_numbers[]);
SHO_PUBLIC void append_units();
SHO_PUBLIC void create_anastomosis(int node_in, int node_out, double anast_radius);
SHO_PUBLIC void define_1d_element_geometry(const char *ELEMFILE);
SHO_PUBLIC void define_1d_element_placenta(const char *ELEMFILE, int anastomosis_elem_in);
SHO_PUBLIC void define_node_geometry(const char *NODEFILE);
SHO_PUBLIC void define_rad_from_file(const char *FIELDFILE, const char *order_system, double s_ratio);
SHO_PUBLIC void define_rad_from_geom(const char *ORDER_SYSTEM, double CONTROL_PARAM, const char *START_FROM,
                                     double START_RAD, const char *GROUP_TYPE, const char *GROUP_OPTIONS);
SHO_PUBLIC void define_ven_rad_from_art(const char *FILENAME, double factor);
SHO_PUBLIC void define_capillary_model(int define_convolutes, int define_generations, int define_parallel,const char *DEFINE_MODEL);
SHO_PUBLIC void element_connectivity_1d();
SHO_PUBLIC void evaluate_ordering();
SHO_PUBLIC void update_radius_by_order(int order,double new_radius,const char *UPDATE_TYPE, double control_parameter);
SHO_PUBLIC void update_1d_elem_field(int ne_field, int elem_number, double value);

#endif /* REPROSIM_GEOMETRY_H */
