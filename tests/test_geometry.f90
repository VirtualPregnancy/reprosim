
module test_geometry
  use testdrive, only : new_unittest, unittest_type, error_type, check
  implicit none
  private

  public :: collect_geometry

contains

!> Collect all exported unit tests
subroutine collect_geometry(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)
  testsuite = [ &
    new_unittest("test_define_node_geometry", test_define_node_geometry), &
    new_unittest("test_define_1d_element_placenta", test_define_1d_element_placenta), &
    new_unittest("test_append_units", test_append_units), &
    new_unittest("test_add_matching_mesh", test_add_matching_mesh), &
    new_unittest("test_define_rad_from_geom", test_define_rad_from_geom) &
    ]

end subroutine collect_geometry

   subroutine test_define_node_geometry(error)
      use arrays,only: dp,node_xyz,num_nodes,nodes
      use geometry,only: define_node_geometry
      use other_consts, only: MAX_FILENAME_LEN
      use indices, only: perfusion_indices
      use test_data, only: write_node_file,delete_node_file
      implicit none
      
      type(error_type), allocatable, intent(out) :: error
      character(len=MAX_FILENAME_LEN) :: NODEFILE
      real(dp) :: test_node_xyz(3,4) ! (3:num_nodes)
      integer :: test_nodes(4) !num_nodes

      NODEFILE = ""
      call write_node_file(NODEFILE)

      !set up indices
      call perfusion_indices()
      !call the subroutine with the test node file
      call define_node_geometry(NODEFILE)
      
      !populate test variables - variables populated by define_node_geometry

      call set_nodes(test_nodes)
      call set_node_xyz(test_node_xyz)  
     
      call check(error, 4, num_nodes)
      if (allocated(error)) return

      call check(error, test_nodes, nodes)
      if (allocated(error)) return

      call check(error, test_node_xyz, node_xyz)
      if (allocated(error)) return
      
      call delete_node_file(NODEFILE) 

   end subroutine test_define_node_geometry

! subroutines to populate test data

  subroutine set_nodes(nodes)
    implicit none 
    integer, intent(out) :: nodes(4)    
    nodes(1)=1
    nodes(2)=2
    nodes(3)=3
    nodes(4)=4    
  end subroutine set_nodes
  
  subroutine set_node_xyz(node_xyz)
    use arrays, only:dp
    implicit none
    real(dp),intent(out) :: node_xyz(3,4)
    node_xyz(1,1) = 0.0_dp
    node_xyz(2,1) = 0.0_dp
    node_xyz(3,1) = 0.0_dp
    node_xyz(1,2) = 0.0_dp 
    node_xyz(2,2) = 0.0_dp
    node_xyz(3,2) = -100.0_dp 
    node_xyz(1,3) = 0.0_dp
    node_xyz(2,3) = -50.0_dp
    node_xyz(3,3) = -150.0_dp
    node_xyz(1,4) = 0.0_dp
    node_xyz(2,4) = 50.0_dp
    node_xyz(3,4) = -150.0_dp  
  end subroutine set_node_xyz

!***************************************************

   subroutine test_define_1d_element_placenta(error)
      use arrays,only: dp,num_elems,elem_nodes,elems_at_node,elem_cnct, &
           num_nodes,node_xyz,nodes,elem_field,elem_direction,elem_ordrs,maxgen
      use geometry,only: define_1d_element_placenta, define_node_geometry
      use indices, only: ne_length,num_ord,perfusion_indices
      use other_consts, only: MAX_FILENAME_LEN
      use test_data, only: write_elem_file,delete_elem_file, write_node_file, delete_node_file
      use diagnostics, only: get_diagnostics_level, set_diagnostics_level
      implicit none

      type(error_type), allocatable, intent(out) :: error
      character(len=MAX_FILENAME_LEN) :: ELEMFILE, NODEFILE
      integer :: test_elem_nodes(2,3) !(2,num_elems)
      integer :: test_elems_at_node(4,0:3)    
      integer :: test_elem_cnct(-1:1,0:2,0:3) !array size: -1:1,0:2,0:num_elems    
      real(dp) :: test_elem_field(1,3) !(1,num_elems) 
      real(dp) :: test_elem_direction(3,3) !array size: 3,num_elems   
      integer :: test_elem_ordrs(num_ord,3) !array size: num_ord, num_elems 

      call set_diagnostics_level(0)
      ELEMFILE = ""
      NODEFILE = ""

      call write_node_file(NODEFILE)
      call write_elem_file(ELEMFILE)
          
      call define_node_geometry(NODEFILE)
      call define_1d_element_placenta(ELEMFILE)
     
      !populate test variables - variables populated by sub define_1d_element_placenta
      call set_elem_nodes(test_elem_nodes)
      call set_elems_at_node(test_elems_at_node)
      call set_elem_cnct(test_elem_cnct)  
      call set_elem_length(test_elem_field)
      call set_elem_direction(test_elem_direction) 
      call set_elem_ordrs(test_elem_ordrs)  
    
      call check(error, 3, num_elems)
      if (allocated(error)) return

      call check(error, test_elem_nodes, elem_nodes)
      if (allocated(error)) return
      
      call check(error, test_elems_at_node, elems_at_node)
      if (allocated(error)) return

      call check(error, test_elem_cnct, elem_cnct)
      if (allocated(error)) return

      call check(error, test_elem_field(1,3), elem_field(ne_length,num_elems))
      if (allocated(error)) return

      call check(error, test_elem_direction, elem_direction)
      if (allocated(error)) return

      call check(error, test_elem_ordrs, elem_ordrs)
      if (allocated(error)) return

      call check(error, 2, maxgen)
      if (allocated(error)) return
      
      call delete_elem_file(ELEMFILE)
      call delete_node_file(NODEFILE)
     
   end subroutine test_define_1d_element_placenta

! subroutines to populate test data 
  
  subroutine set_elem_nodes(elem_nodes)
     implicit none
     integer,intent(out) :: elem_nodes(2,3)  
     elem_nodes(1,1)=1
     elem_nodes(2,1)=2
     elem_nodes(1,2)=2
     elem_nodes(2,2)=3
     elem_nodes(1,3)=2
     elem_nodes(2,3)=4
  end subroutine set_elem_nodes
  
  subroutine set_elems_at_node(elems_at_node)
     implicit none
     integer,intent(out) :: elems_at_node(4,0:3) !array size: num_nodes, 0:3
     elems_at_node=0
     elems_at_node(1,0)=1
     elems_at_node(1,1)=1
     elems_at_node(2,0)=3
     elems_at_node(2,1)=1
     elems_at_node(2,2)=2
     elems_at_node(2,3)=3
     elems_at_node(3,0)=1
     elems_at_node(3,1)=2
     elems_at_node(4,0)=1
     elems_at_node(4,1)=3    
  end subroutine set_elems_at_node
    
  subroutine set_elem_cnct(elem_cnct)
    implicit none
    integer,intent(out) :: elem_cnct(-1:1,0:2,0:3) !array size: -1:1,0:2,0:num_elems
    elem_cnct = 0
    elem_cnct(1,0,1)=2
    elem_cnct(1,1,1)=2
    elem_cnct(1,2,1)=3
    elem_cnct(-1,0,2)=1
    elem_cnct(-1,1,2)=1
    elem_cnct(-1,0,3)=1
    elem_cnct(-1,1,3)=1      
  end subroutine set_elem_cnct    
  
  subroutine set_elem_length(elem_field)
    use arrays, only:dp
    implicit none
    real(dp),intent(inout) :: elem_field(1,3) ! array size: num_ne, num_elems

    elem_field(1,1) = 100.0_dp
    elem_field(1,2) = 70.710678118654755_dp
    elem_field(1,3) = 70.710678118654755_dp  
  end subroutine set_elem_length
  
  subroutine set_elem_direction(elem_direction)
    use arrays, only:dp
    implicit none
    real(dp),intent(out) :: elem_direction(3,3) !array size: 3,num_elems  
    elem_direction(1,1) = 0.0_dp
    elem_direction(2,1) = 0.0_dp
    elem_direction(3,1) = -1.0_dp
    elem_direction(1,2) = 0.0_dp
    elem_direction(2,2)= -0.70710678118654746_dp
    elem_direction(3,2) = -0.70710678118654746_dp
    elem_direction(1,3) = 0.0_dp
    elem_direction(2,3) = 0.70710678118654746_dp
    elem_direction(3,3) = -0.70710678118654746_dp  
  end subroutine set_elem_direction


  subroutine set_elem_ordrs(elem_ordrs)
    use indices, only: num_ord
    implicit none
    integer,intent(out) :: elem_ordrs(num_ord,3) !array size: num_ord, num_elems
    elem_ordrs(1,1) = 1
    elem_ordrs(2,1) = 2
    elem_ordrs(3,1) = 2
    elem_ordrs(1,2) = 2
    elem_ordrs(2,2) = 1
    elem_ordrs(3,2) = 1
    elem_ordrs(1,3) = 2
    elem_ordrs(2,3) = 1
    elem_ordrs(3,3) = 1 
  end subroutine set_elem_ordrs
  
!***************************************************
  
  subroutine test_append_units(error)
    use arrays,only: num_units,units,elem_units_below
    use indices,only: num_nu
    use geometry, only:append_units
    implicit none

    type(error_type), allocatable, intent(out) :: error
    integer :: test_units(2)
    integer :: test_elem_units_below(1:3) ! 1:num_elems
      
    call append_units()
    
    !populate test variables  
    call set_units(test_units)
    call set_elem_units_below(test_elem_units_below)
            
    call check(error, 2, num_units)
    if (allocated(error)) return

    call check(error, test_units, units)
    if (allocated(error)) return

    call check(error, test_elem_units_below, elem_units_below)
    if (allocated(error)) return
            
  end subroutine test_append_units

!subroutines to populate test data 
  subroutine set_units(units)
    implicit none
    integer,intent(out) :: units(2)
    units(1)=2
    units(2)=3
  end subroutine set_units
   
  subroutine set_elem_units_below(elem_units_below)
    implicit none
    integer,intent(out) :: elem_units_below(1:3) ! 1:num_elems
    elem_units_below(1) = 2 
    elem_units_below(2) = 1 
    elem_units_below(3) = 1
  end subroutine set_elem_units_below

!***************************************************
   
   subroutine test_add_matching_mesh(error)
      use arrays,only: dp,num_nodes,num_elems,nodes,node_xyz,elems,elem_nodes,&
                       elems_at_node,elem_field,elem_direction,elem_cnct,elem_ordrs                      
      use geometry, only:add_matching_mesh
      use indices, only: num_ord
      implicit none

      type(error_type), allocatable, intent(out) :: error
      integer :: test_nodes(8) !arterial and venous num_nodes
      real(dp) :: test_node_xyz(3,8) !(3:num_nodes)
      integer :: test_elem_nodes(2,8) !(2:num_elems)
      integer :: test_elems_at_node(8,0:3) !array size: num_nodes, 0:3
      integer :: test_elem_cnct(-1:1,0:2,0:8) !array size: -1:1,0:2,0:num_elems
      integer :: test_elem_ordrs(num_ord,8) !array size: num_ord, num_elems
      character(LEN=100) :: umbilical_elem_option     

      umbilical_elem_option = "same_as_arterial"
      call add_matching_mesh(umbilical_elem_option)

      !populate test variables
      call set_nodes(test_nodes(1:4)) !arterial num_nodes
      call set_nodes_ven(test_nodes(5:8)) !venous num_nodes
      call set_node_xyz(test_node_xyz(1:3,1:4))
      call set_node_xyz_ven(test_node_xyz(1:3,5:8))
      call set_elem_nodes(test_elem_nodes(1:2,1:3))
      call set_elem_nodes_ven(test_elem_nodes(1:2,4:8))
      call set_elems_at_node_new(test_elems_at_node)
      call set_elem_cnct_new(test_elem_cnct)
      call set_elem_ordrs(test_elem_ordrs(1:3,1:3))
      call set_elem_ordrs_ven(test_elem_ordrs(1:3,4:8))
      
      call check(error, 8, num_nodes)
      if (allocated(error)) return

      call check(error, 8, num_elems)
      if (allocated(error)) return

      call check(error, test_nodes, nodes)
      if (allocated(error)) return

      call check(error, test_node_xyz, node_xyz)
      if (allocated(error)) return

      call check(error, test_elem_nodes, elem_nodes)
      if (allocated(error)) return

      call check(error, test_elems_at_node, elems_at_node)
      if (allocated(error)) return

      call check(error, test_elem_cnct, elem_cnct)
      if (allocated(error)) return

      call check(error, test_elem_ordrs, elem_ordrs)
      if (allocated(error)) return
         
   end subroutine test_add_matching_mesh

 !subroutines to populate test data   
   subroutine set_nodes_ven(nodes) 
    implicit none
    integer, intent(out) :: nodes(4)    
    nodes(1)=5
    nodes(2)=6
    nodes(3)=7
    nodes(4)=8    
   end subroutine set_nodes_ven
   
   subroutine set_node_xyz_ven(node_xyz)
    use arrays, only:dp
    implicit none
    real(dp),intent(out) :: node_xyz(3,4)
    node_xyz(1,1) = 0.0_dp
    node_xyz(2,1) = 9.9999999999999995E-007_dp
    node_xyz(3,1) = 0.0_dp
    node_xyz(1,2) = 0.0_dp 
    node_xyz(2,2) = 9.9999999999999995E-007_dp
    node_xyz(3,2) = -100.0_dp 
    node_xyz(1,3) = 0.0_dp
    node_xyz(2,3) = -49.999999000000003_dp
    node_xyz(3,3) = -150.0_dp
    node_xyz(1,4) = 0.0_dp
    node_xyz(2,4) = 50.000000999999997_dp
    node_xyz(3,4) = -150.0_dp  
  end subroutine set_node_xyz_ven
  
  subroutine set_elem_nodes_ven(elem_nodes)
     implicit none
     integer,intent(out) :: elem_nodes(2,5)  
     elem_nodes(1,1)=6
     elem_nodes(2,1)=5
     elem_nodes(1,2)=7
     elem_nodes(2,2)=6
     elem_nodes(1,3)=8
     elem_nodes(2,3)=6
     elem_nodes(1,4)=3
     elem_nodes(2,4)=7
     elem_nodes(1,5)=4
     elem_nodes(2,5)=8
  end subroutine set_elem_nodes_ven
  
  subroutine set_elems_at_node_new(elems_at_node)
     implicit none
     integer,intent(out) :: elems_at_node(8,0:3) !array size: num_nodes, 0:3
     elems_at_node=0
     elems_at_node(1,0)=1
     elems_at_node(1,1)=1
     elems_at_node(2,0)=3
     elems_at_node(2,1)=1
     elems_at_node(2,2)=2
     elems_at_node(2,3)=3
     elems_at_node(3,0)=2
     elems_at_node(3,1)=2
     elems_at_node(3,2)=7
     elems_at_node(4,0)=2
     elems_at_node(4,1)=3  
     elems_at_node(4,2)=8  
     elems_at_node(5,0)=1
     elems_at_node(5,1)=4     
     elems_at_node(6,0)=3
     elems_at_node(6,1)=4
     elems_at_node(6,2)=5
     elems_at_node(6,3)=6   
     elems_at_node(7,0)=2
     elems_at_node(7,1)=5  
     elems_at_node(7,2)=7
     elems_at_node(8,0)=2
     elems_at_node(8,1)=6  
     elems_at_node(8,2)=8      
  end subroutine set_elems_at_node_new
    
  subroutine set_elem_cnct_new(elem_cnct)
    implicit none
    integer,intent(out) :: elem_cnct(-1:1,0:2,0:8) !array size: -1:1,0:2,0:num_elems
    elem_cnct = 0
    elem_cnct(1,0,1)=2
    elem_cnct(1,1,1)=2
    elem_cnct(1,2,1)=3
    elem_cnct(-1,0,2)=1
    elem_cnct(-1,1,2)=1
    elem_cnct(1,0,2)=1
    elem_cnct(1,1,2)=7   
    elem_cnct(-1,0,3)=1
    elem_cnct(-1,1,3)=1    
    elem_cnct(1,0,3)=1
    elem_cnct(1,1,3)=8  
    elem_cnct(-1,0,4)=2
    elem_cnct(-1,1,4)=5 
    elem_cnct(-1,2,4)=6     
    elem_cnct(-1,0,5)=1
    elem_cnct(-1,1,5)=7  
    elem_cnct(1,0,5)=1
    elem_cnct(1,1,5)=4  
    elem_cnct(-1,0,6)=1
    elem_cnct(-1,1,6)=8 
    elem_cnct(1,0,6)=1
    elem_cnct(1,1,6)=4
    elem_cnct(-1,0,7)=1
    elem_cnct(-1,1,7)=2
    elem_cnct(1,0,7)=1
    elem_cnct(1,1,7)=5
    elem_cnct(-1,0,8)=1
    elem_cnct(-1,1,8)=3
    elem_cnct(1,0,8)=1
    elem_cnct(1,1,8)=6
  end subroutine set_elem_cnct_new    
  
  subroutine set_elem_ordrs_ven(elem_ordrs)
    use indices, only: num_ord
    implicit none
    integer,intent(out) :: elem_ordrs(num_ord,5)
    elem_ordrs(1,1) = 1
    elem_ordrs(2,1) = 2
    elem_ordrs(3,1) = 2
    elem_ordrs(1,2) = 2
    elem_ordrs(2,2) = 1
    elem_ordrs(3,2) = 1   
    elem_ordrs(1,3) = 2
    elem_ordrs(2,3) = 1
    elem_ordrs(3,3) = 1     
    elem_ordrs(1,4) = 2
    elem_ordrs(2,4) = 1
    elem_ordrs(3,4) = 1 
    elem_ordrs(1,5) = 2
    elem_ordrs(2,5) = 1
    elem_ordrs(3,5) = 1
  end subroutine set_elem_ordrs_ven
  
!***************************************************
  
  subroutine test_define_rad_from_geom(error)
    use arrays,only: dp,elem_field,num_elems                     
    use geometry, only:define_rad_from_geom
    use indices, only: ne_radius
    implicit none  
  
    type(error_type), allocatable, intent(out) :: error
    character(LEN=100) :: order_system,order_options,name
    real(dp) :: s_ratio, inlet_rad
    real(dp) :: test_elem_field(1,8) !(1: num_elems)

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
  
    call set_elem_radius(test_elem_field)

    
    call check(error, test_elem_field(1,1:8), elem_field(ne_radius,1:num_elems))
    if (allocated(error)) return

 
  end subroutine test_define_rad_from_geom

 !subroutines to populate test data
  subroutine set_elem_radius(elem_f)
    use arrays, only: dp
    implicit none
    real(dp),intent(out) :: elem_f(1,8) !(1 : num_elems)
    
    elem_f(1,1)=3.0_dp
    elem_f(1,2)=1.9480519480519483_dp
    elem_f(1,3)=1.9480519480519483_dp
    elem_f(1,4)=5.0000000000000009_dp
    elem_f(1,5)=3.2258064516129035_dp
    elem_f(1,6)=3.2258064516129035_dp
    elem_f(1,7)=3.2258064516129035_dp
    elem_f(1,8)=3.2258064516129035_dp
  end subroutine set_elem_radius
  
end module test_geometry

