module test_prf_inlet_flow
  use testdrive, only : new_unittest, unittest_type, error_type, check
  implicit none
  private

  public :: collect_prf_inlet_flow

contains

!> Collect all exported unit tests
subroutine collect_prf_inlet_flow(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)
  testsuite = [ &
    new_unittest("test_flow", test_flow) &
    ]

end subroutine collect_prf_inlet_flow

  subroutine test_flow(error)
    use arrays, only: dp,elem_field,node_field,num_elems,num_nodes, &
                      num_units,unit_field 
    use indices, only: ne_Qdot,nj_bv_press,nu_perf,nu_blood_press,perfusion_indices
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: define_node_geometry,define_1d_element_placenta,append_units, &
                        add_matching_mesh,define_rad_from_geom,define_capillary_model
    use pressure_resistance_flow, only: evaluate_prq
    use test_data
    implicit none

    type(error_type), allocatable, intent(out) :: error
    character(len=MAX_FILENAME_LEN) :: NODEFILE
    character(len=MAX_FILENAME_LEN) :: ELEMFILE
    character(LEN=100) :: order_system,order_options,name
    real(dp) :: s_ratio, inlet_rad
    character(LEN=100) :: mesh_type, bc_type, rheology_type,vessel_type,umbilical_elem_option
    real(dp) :: inlet_flow,inlet_pressure,outlet_pressure
    real(dp) :: test_elem_field(1,8) !(num_ne: num_elems)
    real(dp) :: test_node_field(1,8) !(num_nj : num_nodes)
    real(dp) :: test_unit_field(2,2) !(2 : num_units)
    character(LEN=60) :: cap_model
   
    !call subroutines which need to be executed before evaluate_prq   
    NODEFILE = ""
    call write_node_file(NODEFILE)    
    call perfusion_indices()
    call define_node_geometry(NODEFILE) 
    call delete_node_file(NODEFILE) 
    ELEMFILE = ""
    call write_elem_file(ELEMFILE)          
    call define_1d_element_placenta(ELEMFILE)
    call delete_elem_file(ELEMFILE)
    call append_units()
    umbilical_elem_option = "same_as_arterial"
    call add_matching_mesh(umbilical_elem_option)

    !define radius for arterial vessels
    order_system = "strahler"
    s_ratio=1.54_dp
    name = "inlet"
    inlet_rad=3.0_dp    
    order_options = "arterial"   
    call define_rad_from_geom(order_system, s_ratio, name, inlet_rad, order_options)
    !define radius for venous vessels
    order_system = "strahler"
    s_ratio=1.55_dp
    name = ""
    inlet_rad=5.0_dp  
    order_options = "venous"
    call define_rad_from_geom(order_system, s_ratio, name, inlet_rad, order_options)
    cap_model = 'byrne_simplified'
    call define_capillary_model(10,3,6,cap_model)
      
    mesh_type = "full_plus_tube"
    bc_type = "flow"
    rheology_type = "constant_visc"
    vessel_type = "rigid"
    inlet_flow = 111666.7_dp
    inlet_pressure = 0.0_dp
    outlet_pressure = 2660.0_dp
    call evaluate_prq(mesh_type,bc_type,rheology_type,vessel_type,inlet_flow,inlet_pressure,outlet_pressure)
      
    !populate test data
    call set_flow(test_elem_field)
    
    !check flow
    call check(error, test_elem_field(1,1:8),elem_field(ne_Qdot,1:num_elems))
    if (allocated(error)) return

  end subroutine test_flow

!subroutines to populate test data

  subroutine set_flow(elem_field)
    use arrays, only: dp
    implicit none
    real(dp),intent(out) :: elem_field(1,8) !(num_ne : num_elems)
    
    elem_field(1,1) = 111666.7_dp 
    elem_field(1,2) = 55833.350000000020_dp
    elem_field(1,3) = 55833.349999999969_dp
    elem_field(1,4) = 111666.69999999992_dp
    elem_field(1,5) = 55833.350000000006_dp
    elem_field(1,6) = 55833.349999999969_dp
    elem_field(1,7) = 55833.350000000006_dp
    elem_field(1,8) = 55833.349999999969_dp
    
  end subroutine set_flow

                                  
end module test_prf_inlet_flow
