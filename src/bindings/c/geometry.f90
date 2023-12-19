module geometry_c
  implicit none
  private

contains
!
!###################################################################################
!
!*add_matching_mesh:* Replicates an existing mesh, continuing node and element numbers
  subroutine add_matching_mesh_c(umbilical_elem_option,umbilical_elem_option_len, &
        umbilical_element_numbers, umbilical_element_numbers_len) bind(C, name="add_matching_mesh_c")
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_STRING_LEN,MAX_FILENAME_LEN
    use geometry, only: add_matching_mesh
    implicit none

    integer,intent(in) :: umbilical_elem_option_len, umbilical_element_numbers_len
    type(c_ptr), value, intent(in) :: umbilical_elem_option
    character(len=MAX_STRING_LEN) :: umbilical_elem_option_f
    integer,intent(in) :: umbilical_element_numbers(umbilical_element_numbers_len)

    call strncpy(umbilical_elem_option_f, umbilical_elem_option, umbilical_elem_option_len)
#if defined _WIN32 && defined __INTEL_COMPILER
    call so_add_matching_mesh(umbilical_elem_option_f, umbilical_element_numbers)
#else
    call add_matching_mesh(umbilical_elem_option_f, umbilical_element_numbers)
#endif

  end subroutine add_matching_mesh_c

!
!###################################################################################
!
!*append_units:* Appends terminal units at the end of a tree structure
  subroutine append_units_c() bind(C, name="append_units_c")
    use geometry, only: append_units
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_append_units
#else
    call append_units
#endif

  end subroutine append_units_c


!
!###################################################################################
!
!
  subroutine create_anastomosis_c(node_in, node_out, anast_radius) &
    bind(C, name="create_anastomosis_c")
    use geometry, only: create_anastomosis
    use arrays, only: dp

    integer, intent(in) :: node_in,node_out
    real(dp), intent(in) :: anast_radius

#if defined _WIN32 && defined __INTEL_COMPILER
   call so_create_anastomosis(node_in, node_out, anast_radius)
#else
  call create_anastomosis(node_in, node_out, anast_radius)
#endif

  end subroutine create_anastomosis_c

!
!###################################################################################
!
  subroutine define_1d_element_geometry_c(ELEMFILE, filename_len) &
              bind(C, name="define_1d_element_geometry_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: define_1d_element_geometry
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: ELEMFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, ELEMFILE, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_1d_element_geometry(filename_f)
#else
    call define_1d_element_geometry(filename_f)
#endif

  end subroutine define_1d_element_geometry_c

!
!###################################################################################
!
  subroutine define_1d_element_placenta_c(ELEMFILE, filename_len, anastomosis_elem_in) &
              bind(C, name="define_1d_element_placenta_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: define_1d_element_placenta
    implicit none

    integer,intent(in) :: filename_len
    integer,intent(in) :: anastomosis_elem_in
    type(c_ptr), value, intent(in) :: ELEMFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, ELEMFILE, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_1d_element_placenta(filename_f,anastomosis_elem_in)
#else
    call define_1d_element_placenta(filename_f,anastomosis_elem_in)
#endif

  end subroutine define_1d_element_placenta_c
!
!###################################################################################
!
  subroutine define_node_geometry_c(NODEFILE, filename_len) bind(C, name="define_node_geometry_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: define_node_geometry
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: NODEFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, NODEFILE, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_node_geometry(filename_f)
#else
    call define_node_geometry(filename_f)
#endif

  end subroutine define_node_geometry_c
!
!###################################################################################
!
  subroutine define_rad_from_file_c(FIELDFILE, filename_len, order_system, &
          order_system_len,s_ratio) bind(C, name="define_rad_from_file_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
    use geometry, only: define_rad_from_file
    use arrays, only: dp
    implicit none

    integer,intent(in) :: filename_len, order_system_len
    type(c_ptr), value, intent(in) :: FIELDFILE, order_system
    real(dp),intent(in) :: s_ratio
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: order_system_f

    call strncpy(filename_f, FIELDFILE, filename_len)
    call strncpy(order_system_f, order_system, order_system_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_rad_from_file(filename_f, order_system_f,s_ratio)
#else
    call define_rad_from_file(filename_f, order_system_f,s_ratio)
#endif

    end subroutine define_rad_from_file_c
!
!##################################################################################
!
!*define_rad_from_geom:* Defines vessel radius based on their geometric structure
  subroutine define_rad_from_geom_c(order_system, order_system_len, control_param, &
        start_from, start_from_len, start_rad, group_type, group_type_len, group_options, group_options_len) &
        bind(C, name="define_rad_from_geom_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_STRING_LEN
    use arrays, only: dp
    use geometry, only: define_rad_from_geom
    implicit none

    real(dp),intent(in) :: control_param, start_rad
    integer,intent(in) :: order_system_len, start_from_len, group_type_len, group_options_len
    type(c_ptr), value, intent(in) :: order_system, start_from, group_type, group_options
    character(len=MAX_STRING_LEN) :: order_system_f, start_from_f, group_type_f, group_options_f

    call strncpy(order_system_f, order_system, order_system_len)
    call strncpy(start_from_f, start_from, start_from_len)
    call strncpy(group_options_f, group_options, group_options_len)
    call strncpy(group_type_f, group_type, group_type_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_rad_from_geom(order_system_f, control_param, start_from_f, start_rad, group_type_f, group_options_f)
#else
    call define_rad_from_geom(order_system_f, control_param, start_from_f, start_rad, group_type_f, group_options_f)
#endif

  end subroutine define_rad_from_geom_c
!
!###########################################################################
!
  subroutine define_ven_rad_from_art_c(FILENAME, filename_len,&
          factor) bind(C, name="define_ven_rad_from_art_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: define_ven_rad_from_art
    use arrays, only: dp
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: FILENAME
    real(dp),intent(in) :: factor
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, FILENAME, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_ven_rad_from_art(filename_f, factor)
#else
    call define_ven_rad_from_art(filename_f,factor)
#endif

    end subroutine define_ven_rad_from_art_c

!
!#################################################################################
!
 subroutine define_capillary_model_c(define_convolutes,define_generations,define_parallel,&
          define_model,define_model_len) bind(C, name="define_capillary_model_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: define_capillary_model
    use arrays, only: dp
    implicit none
    integer, intent(in) :: define_convolutes,define_generations,define_parallel
    integer,intent(in) :: define_model_len
    type(c_ptr), value, intent(in) :: define_model
    character(len=MAX_FILENAME_LEN) :: define_model_f

    call strncpy(define_model_f, define_model, define_model_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_capillary_model(define_convolutes,define_generations,define_parallel,&
          define_model_f)
#else
    call define_capillary_model(define_convolutes,define_generations,define_parallel,&
          define_model_f)
#endif

    end subroutine define_capillary_model_c
!
!##################################################################################
!
!*element_connectivity_1d:*  Calculates element connectivity in 1D and stores in elelem_cnct
  subroutine element_connectivity_1d_c() bind(C, name="element_connectivity_1d_c")
    use geometry, only: element_connectivity_1d
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_element_connectivity_1d
#else
    call element_connectivity_1d
#endif

  end subroutine element_connectivity_1d_c

!
!###################################################################################
!
!*evaluate_ordering:* calculates generations, Horsfield orders, Strahler orders for a given tree
  subroutine evaluate_ordering_c() bind(C, name="evaluate_ordering_c")
    use geometry, only: evaluate_ordering
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_evaluate_ordering
#else
    call evaluate_ordering
#endif

  end subroutine evaluate_ordering_c


!
!###################################################################################
!
subroutine update_1d_elem_field_c(ne_field, elem_number,value) &
   bind(C, name="update_1d_elem_field_c")
    use geometry, only: update_1d_elem_field
    use arrays, only: dp
    implicit none
    integer, intent(in) :: ne_field
    integer, intent(in) :: elem_number
    real(dp), intent(in) :: value

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_update_1d_elem_field(ne_field,elem_number,value)
#else
    call update_1d_elem_field(ne_field,elem_number,value)
#endif
end subroutine update_1d_elem_field_c

end module geometry_c

