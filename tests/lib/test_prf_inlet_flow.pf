@Test
  subroutine test_prf_inlet_flow()
    use pfunit_mod
    use arrays, only: dp,elem_field,node_field,num_elems,num_nodes, &
                      num_units,unit_field 
    use indices, only: ne_Qdot,nj_bv_press,nu_perf,nu_blood_press,perfusion_indices
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: define_node_geometry,define_1d_elements,append_units, &
                        add_matching_mesh,define_rad_from_geom
    use pressure_resistance_flow, only: evaluate_prq
    use test_data
    implicit none

    character(len=MAX_FILENAME_LEN) :: NODEFILE
    character(len=MAX_FILENAME_LEN) :: ELEMFILE
    character(LEN=100) :: order_system,order_options,name
    real(dp) :: s_ratio, inlet_rad
    character(LEN=60) :: mesh_type, bc_type, rheology_type,vessel_type,umbilical_elem_option
    real(dp) :: inlet_flow,inlet_pressure,outlet_pressure
    real(dp) :: test_elem_field(1,8) !(num_ne: num_elems)
    real(dp) :: test_node_field(1,8) !(num_nj : num_nodes)
    real(dp) :: test_unit_field(2,2) !(2 : num_units)
   
    !call subroutines which need to be executed before evaluate_prq   
    NODEFILE = ""
    call write_node_file(NODEFILE)    
    call perfusion_indices()
    call define_node_geometry(NODEFILE) 
    call delete_node_file(NODEFILE) 
    ELEMFILE = ""
    call write_elem_file(ELEMFILE)          
    call define_1d_elements(ELEMFILE)
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
    call define_rad_from_geom(order_system, s_ratio, name, inlet_rad, order_options,"")   
    !define radius for venous vessels
    order_system = "strahler"
    s_ratio=1.55_dp
    name = ""
    inlet_rad=5.0_dp  
    order_options = "venous"
    call define_rad_from_geom(order_system, s_ratio, name, &
                    inlet_rad, order_options,"")
      
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
    call set_pressure(test_node_field)
    call set_terminal_solution(test_unit_field)
       
    !check flow
    @assertEqual(test_elem_field(1,1:8),elem_field(ne_Qdot,1:num_elems)) 
    !check pressure
    @assertEqual(test_node_field(1,1:8),node_field(nj_bv_press,1:num_nodes)) 
    !check flow and pressure for terminal units 
    @assertEqual(test_unit_field(1,1:2),unit_field(nu_perf,1:num_units)) 
    @assertEqual(test_unit_field(2,1:2),unit_field(nu_blood_press,1:num_units))   

  end subroutine test_prf_inlet_flow

!subroutines to populate test data

  subroutine set_flow(elem_field)
    use arrays, only: dp
    implicit none
    real(dp),intent(out) :: elem_field(1,8) !(num_ne : num_elems)
    
    elem_field(1,1) = 111666.7_dp 
    elem_field(1,2) = 55833.294838280883_dp     
    elem_field(1,3) = 55833.306428283686_dp     
    elem_field(1,4) = 111666.80201897584_dp     
    elem_field(1,5) = 55833.393571716500_dp     
    elem_field(1,6) = 55833.306428283686_dp     
    elem_field(1,7) = 55833.393571716500_dp     
    elem_field(1,8) = 55833.306428283686_dp 
    
  end subroutine set_flow

  subroutine set_pressure(node_field)
    use arrays, only: dp
    implicit none
    real(dp),intent(out) :: node_field(1,8) !(1 : num_nodes)
  
    node_field(1,1) = 6650.0037645709381_dp     
    node_field(1,2) = 5470.4463792225870_dp     
    node_field(1,3) = 3124.8349354185775_dp     
    node_field(1,4) = 3124.8190499593911_dp     
    node_field(1,5) = 2660.0000000000000_dp     
    node_field(1,6) = 2812.8706292625357_dp     
    node_field(1,7) = 3124.8349310067424_dp     
    node_field(1,8) = 3124.8339425717932_dp 
                               
  end subroutine set_pressure   
   
  subroutine set_terminal_solution(unit_field)
    use arrays, only: dp
    implicit none
    real(dp),intent(out) :: unit_field(2,2) !(2 : num_units)
    !flow
    unit_field(1,1) = 55833.294838280883_dp  
    unit_field(1,2) = 55833.306428283686_dp 
    !pressure    
    unit_field(2,1) = 3124.8349354185775_dp   
    unit_field(2,2) = 3124.8190499593911_dp

  end subroutine set_terminal_solution            
                                  
