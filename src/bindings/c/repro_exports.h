
#ifndef REPROSIM_REPRO_EXPORTS_H
#define REPROSIM_REPRO_EXPORTS_H

#include "symbol_export.h"

SHO_PUBLIC void export_1d_elem_field(int ne_field, const char *EXELEMFILE, const char *group_name, const char *field_name );
SHO_PUBLIC void export_1d_elem_geometry(const char *EXELEMFILE, const char *name);
SHO_PUBLIC void export_node_field(int nj_field, const char *EXNODEFIELD, const char *name, const char *field_name);
SHO_PUBLIC void export_terminal_perfusion(const char *EXNODEFILE, const char *name);
SHO_PUBLIC void export_node_geometry(const char *EXNODEFILE, const char *name);

#endif /* REPROSIM_REPRO_EXPORTS_H */
