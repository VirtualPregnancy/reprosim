#include "pressure_resistance_flow.h"
#include "string.h"

void calculate_stats_c(const char *FLOW_GEN_FILE, int *filename_len, double *image_voxel_size);
void evaluate_prq_c(const char *mesh_type, int *mesh_type_len, const char *bc_type, int *bc_type_len, double *inlet_flow, double *inlet_pressure, double *outlet_pressure);


void calculate_stats(const char *FLOW_GEN_FILE, double image_voxel_size)
{
   int filename_len = strlen(FLOW_GEN_FILE);
   calculate_stats_c(FLOW_GEN_FILE, &filename_len, &image_voxel_size);
}

void evaluate_prq(const char *mesh_type, const char *bc_type, double inlet_flow,double inlet_pressure, double outlet_pressure)
{
	
	int mesh_type_len = strlen(mesh_type);
	int bc_type_len = strlen(bc_type);
    evaluate_prq_c(mesh_type, &mesh_type_len, bc_type, &bc_type_len, &inlet_flow, &inlet_pressure, &outlet_pressure);
}
