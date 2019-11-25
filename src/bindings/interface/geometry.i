
%module(package="reprosim") geometry
%include symbol_export.h

%typemap(in) (int umbilical_element_numbers_len, int umbilical_element_numbers[]) {
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    SWIG_fail;
  }
  $1 = PyList_Size($input);
  $2 = (int *) malloc(($1)*sizeof(int));
  for (i = 0; i < $1; i++) {
    PyObject *o = PyList_GetItem($input, i);
    if (!PyInt_Check(o)) {
      free($2);
      PyErr_SetString(PyExc_ValueError, "List items must be integers");
      SWIG_fail;
    }
    $2[i] = PyInt_AsLong(o);
  }
}

%typemap(freearg) (int umbilical_element_numbers_len, int umbilical_element_numbers[]) {
  if ($2) free($2);
}

%{
#include "geometry.h"
%}

// define_rad_from_file has an optional argument that C cannot replicate,
// so we use SWIG to override with a C++ version that can.
void define_1d_elements(const char *ELEMFILE, int anastomosis_elem_in = 0);
void define_rad_from_file(const char *FIELDFILE, const char *order_system="strahler", double s_ratio=1.54);
void define_rad_from_geom(const char *ORDER_SYSTEM, double CONTROL_PARAM, const char *START_FROM, double START_RAD, const char *group_type_in="all", const char *group_option_in="dummy");
void define_ven_rad_from_art(const char *FILENAME, double factor=2.00);

%include geometry.h

