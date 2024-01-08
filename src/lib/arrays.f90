module arrays
!*Description:* This module defines arrays.
!

  implicit none

  integer :: num_elems,num_nodes,num_units,maxgen,num_arterial_elems, &
             num_conv,num_conv_gen,anastomosis_elem,num_inlets,num_outlets

  integer :: num_elems_fetal,num_nodes_fetal
  
  integer, parameter :: dp=kind(0.d0) !  for double precision
  real(dp),parameter :: zero_tol = 1.0e-12_dp
  real(dp),parameter :: loose_tol = 1.0e-6_dp

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
  integer, allocatable :: umbilical_inlets(:) !arterial inlets
  integer, allocatable :: umbilical_outlets(:) !venous outlet
  integer :: min_art,max_art,min_ven,max_ven

  integer,allocatable :: nodes_fetal(:) !allocated in define_node_geometry
  integer,allocatable :: elems_fetal(:) !allocated in define_1d_elements
  integer,allocatable :: elem_cnct_fetal(:,:,:)  !NXI(-ni:ni,1,ne)
  integer,allocatable :: elem_nodes_fetal(:,:)
  integer,allocatable :: elems_at_node_fetal(:,:)

  real(dp),allocatable :: elem_field(:,:) !properties of elements
  real(dp),allocatable :: elem_direction(:,:)
  real(dp),allocatable :: node_xyz(:,:)
  real(dp),allocatable :: unit_field(:,:) !properties of elastic units
  real(dp),allocatable :: node_field(:,:)

  real(dp),allocatable :: elem_field_fetal(:,:) !properties of elements
  real(dp),allocatable :: elem_direction_fetal(:,:)
  real(dp),allocatable :: node_xyz_fetal(:,:)
  real(dp),allocatable :: node_field_fetal(:,:)

  real(dp) :: cap_resistance,terminal_resistance,terminal_length, &
              cap_radius,total_cap_volume,total_cap_surface_area

  integer :: num_convolutes,num_generations,num_parallel,model_type

  integer :: capillary_model_type

! temporary, for debugging:
  real(dp) :: unit_before

  private
  public set_node_field_value, elem_field, num_elems, elem_nodes, node_xyz, nodes, elems, &
    num_nodes, units, num_units, unit_field, node_field, dp, elem_cnct, elem_ordrs, elem_direction, &
    elems_at_node, elem_symmetry, elem_units_below, maxgen, num_arterial_elems, &
    num_conv,num_conv_gen,cap_resistance,terminal_resistance,terminal_length, &
    cap_radius,elem_cnct_no_anast,anastomosis_elem, &
    is_capillary_unit,total_cap_volume,total_cap_surface_area,umbilical_inlets,umbilical_outlets, &
    art_ven_elem_map,num_inlets,num_outlets,min_art,max_art,min_ven,max_ven,num_convolutes,&
    num_generations,num_parallel,model_type,capillary_model_type,zero_tol,loose_tol

  public elem_field_fetal, num_elems_fetal, elem_nodes_fetal, node_xyz_fetal, nodes_fetal, elems_fetal,&
          num_nodes_fetal,node_field_fetal, elem_cnct_fetal,elem_direction_fetal,elems_at_node_fetal
  
!dec$ attributes dllexport :: elem_nodes, elem_units_below, maxgen, elem_direction, num_elems, num_nodes, elems_at_node
!dec$ attributes dllexport :: elem_ordrs, elem_field, elem_cnct, node_xyz, nodes, units, unit_field, node_field, num_units

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
