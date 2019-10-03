module arrays
!*Description:* This module defines arrays.
!

  implicit none

  integer :: num_elems,num_nodes,num_units,maxgen,num_arterial_elems, &
             num_conv,num_conv_gen,anastomosis_elem
  
  integer, parameter :: dp=kind(0.d0) !  for double precision

  integer,allocatable :: nodes(:) !allocated in define_node_geometry
  integer,allocatable :: elems(:) !allocated in define_1d_elements
  integer,allocatable :: elem_cnct(:,:,:)  !NXI(-ni:ni,1,ne)
  integer,allocatable :: elem_cnct_no_anast(:,:,:)  !NXI(-ni:ni,1,ne) !copy of elem_cnct array without the anastomosis element
  integer,allocatable :: elem_nodes(:,:)
  integer,allocatable :: elem_ordrs(:,:)
  integer,allocatable :: elem_symmetry(:)
  integer,allocatable :: elem_units_below(:)
  integer,allocatable :: elems_at_node(:,:)
  integer,allocatable :: units(:)
  integer,allocatable :: is_capillary_unit(:)
  integer,allocatable :: art_ven_elem_map(:)
  integer :: umbilical_inlets(2) !arterial inlets
  integer :: umbilical_outlets(2) !venous outlet

  real(dp),allocatable :: elem_field(:,:) !properties of elements
  real(dp),allocatable :: elem_direction(:,:)
  real(dp),allocatable :: node_xyz(:,:)
  real(dp),allocatable :: unit_field(:,:) !properties of elastic units
  real(dp),allocatable :: node_field(:,:)

  real(dp) :: cap_resistance,terminal_resistance,terminal_length, &
              cap_radius,total_cap_volume,total_cap_surface_area

! temporary, for debugging:
  real(dp) :: unit_before

  private
  public set_node_field_value, elem_field, num_elems, elem_nodes, node_xyz, nodes, elems, &
    num_nodes, units, num_units, unit_field, node_field, dp, elem_cnct, elem_ordrs, elem_direction, &
    elems_at_node, elem_symmetry, elem_units_below, maxgen, num_arterial_elems, &
    num_conv,num_conv_gen,cap_resistance,terminal_resistance,terminal_length, &
    cap_radius,elem_cnct_no_anast,anastomosis_elem, &
    is_capillary_unit,total_cap_volume,total_cap_surface_area,umbilical_inlets,umbilical_outlets, &
    art_ven_elem_map

contains
  subroutine set_node_field_value(row, col, value)  
  !*Description:* This subroutine sets the value of a node field
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SET_NODE_FIELD_VALUE" :: SET_NODE_FIELD_VALUE
    integer, intent(in) :: row, col
    real(dp), intent(in) :: value
    

    node_field(row, col) = value
	
  end subroutine set_node_field_value


end module arrays
