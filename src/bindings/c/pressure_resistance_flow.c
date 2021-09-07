#include "pressure_resistance_flow.h"
#include "string.h"

void calculate_stats_c(const char *FLOW_GEN_FILE, int *filename_len, double *image_voxel_size, int *output_level);
void evaluate_prq_c(const char *mesh_type, int *mesh_type_len, const char *bc_type, int *bc_type_len, const char *rheology_type, int *rheology_type_len, const char *vessel_type, int *vessel_type_len, double *inlet_flow,double *inlet_pressure, double *outlet_pressure);

void calculate_stats(const char *FLOW_GEN_FILE, double image_voxel_size, int output_level)
{
   int filename_len = strlen(FLOW_GEN_FILE);
   calculate_stats_c(FLOW_GEN_FILE, &filename_len, &image_voxel_size,&output_level);
}

void evaluate_prq(const char *mesh_type, const char *bc_type, const char *rheology_type, const char *vessel_type, double inlet_flow,double inlet_pressure, double outlet_pressure)
{
	int mesh_type_len = strlen(mesh_type);
	int bc_type_len = strlen(bc_type);
	int rheology_type_len = strlen(rheology_type);
	int vessel_type_len = strlen(vessel_type);
    evaluate_prq_c(mesh_type, &mesh_type_len, bc_type, &bc_type_len, rheology_type, &rheology_type_len, vessel_type, &vessel_type_len, &inlet_flow, &inlet_pressure, &outlet_pressure);
}
