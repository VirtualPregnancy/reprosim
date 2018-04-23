module test_data
!*Description:* This module creates the data required by unit tests. 
! It creates files and populates variables that are used by each tested subroutine
!
!*Contributor(s):* Monika Byrne

  implicit none

  private
  public write_node_file, delete_node_file, get_num_nodes, get_dim_nodes, get_nodes, &
         get_dim_node_xyz, get_node_xyz, write_elem_file, delete_elem_file, &
         get_num_elems, get_dim_elem_nodes,get_elem_nodes,get_dim_elems_at_node, &
         get_elems_at_node,get_dim_elem_cnct,get_elem_cnct,get_dim_elem_field, &
         get_elem_field,get_dim_elem_direction,get_elem_direction, &
         get_dim_elem_ordrs,get_elem_ordrs,get_maxgen,get_num_units, &
         get_dim_units,get_units,get_dim_elem_units_below,get_elem_units_below
           
contains
  subroutine write_node_file(NODEFILE)  
    use other_consts, only: MAX_FILENAME_LEN
    character(len=MAX_FILENAME_LEN), intent(out) :: NODEFILE
    NODEFILE = "test.ipnode"
    open(10, file=NODEFILE, status="replace")
    write(10,*) 'CMISS Version 2.1  ipnode File Version 2'
    write(10,*) 'Heading:'
    write(10,*) 'The number of nodes is [ 61742]:  4'
    write(10,*) 'Number of coordinates [3]: 3'
    write(10,*) 'Do you want prompting for different versions of nj=1 [N]? Y'
    write(10,*) 'Do you want prompting for different versions of nj=2 [N]? Y'
    write(10,*) 'Do you want prompting for different versions of nj=3 [N]? Y'
    write(10,*) 'The number of derivatives for coordinate 1 is [0]: 0'
    write(10,*) 'The number of derivatives for coordinate 2 is [0]: 0'
    write(10,*) 'The number of derivatives for coordinate 3 is [0]: 0'
    write(10,*) 'Node number [  1001]:   1'
    write(10,*) 'The number of versions for nj=1 is [1]:  1'
    write(10,*) 'The Xj(1) coordinate is [ 0.00000E+00]:   0.00'
    write(10,*) 'The number of versions for nj=2 is [1]:  1'
    write(10,*) 'The Xj(2) coordinate is [ 0.00000E+00]:   0.00'
    write(10,*) 'The number of versions for nj=3 is [1]:  1'
    write(10,*) 'The Xj(3) coordinate is [ 0.00000E+00]:   0.00'
    write(10,*) ''
    write(10,*) 'Node number [  1002]:   2'
    write(10,*) 'The number of versions for nj=1 is [1]:  1'
    write(10,*) 'The Xj(1) coordinate is [-0.15641E+03]:   0.00'
    write(10,*) 'The number of versions for nj=2 is [1]:  1'
    write(10,*) 'The Xj(2) coordinate is [-0.10446E+03]:   0.00'
    write(10,*) 'The number of versions for nj=3 is [1]:  1'
    write(10,*) 'The Xj(3) coordinate is [-0.12000E+01]:   -100.00'
    write(10,*) ''
    write(10,*) 'Node number [  1003]:   3'
    write(10,*) 'The number of versions for nj=1 is [1]:  1'
    write(10,*) 'The Xj(1) coordinate is [-0.15641E+03]:   0.00'
    write(10,*) 'The number of versions for nj=2 is [1]:  1'
    write(10,*) 'The Xj(2) coordinate is [-0.14300E+03]:   -50.00'
    write(10,*) 'The number of versions for nj=3 is [1]:  1'
    write(10,*) 'The Xj(3) coordinate is [-0.76198E+02]:   -150.00'
    write(10,*) ''
    write(10,*) 'Node number [  1004]:   4'
    write(10,*) 'The number of versions for nj=1 is [1]:  1'
    write(10,*) 'The Xj(1) coordinate is [-0.18992E+03]:   0.00'
    write(10,*) 'The number of versions for nj=2 is [1]:  1'
    write(10,*) 'The Xj(2) coordinate is [-0.15305E+03]:   50.00'
    write(10,*) 'The number of versions for nj=3 is [1]:  1'
    write(10,*) 'The Xj(3) coordinate is [-0.11340E+03]:   -150.00'
    close(10)	
  end subroutine write_node_file

  subroutine delete_node_file(NODEFILE)
    use other_consts, only: MAX_FILENAME_LEN 
    character(len=MAX_FILENAME_LEN), intent(out) :: NODEFILE  
    NODEFILE = "test.ipnode"
    open(10,file=NODEFILE, status="old")
    close(10, status="delete")     
  end subroutine delete_node_file
  
  subroutine get_num_nodes(tree_type,num_nodes)
    character(len=20),intent(in) :: tree_type
    integer,intent(out) :: num_nodes
    if(tree_type.EQ."arterial")then
      num_nodes = 4
    elseif(tree_type.EQ."arterial_and_venous")then
      num_nodes = 8
    endif
  end subroutine get_num_nodes  
  
  subroutine get_dim_nodes(x)
    integer,intent(out) :: x
    x = 4
  end subroutine get_dim_nodes
  
  subroutine get_nodes(nodes) 
    integer, intent(out) :: nodes(4)    
    nodes(1)=1
    nodes(2)=2
    nodes(3)=3
    nodes(4)=4    
  end subroutine get_nodes
  
  subroutine get_dim_node_xyz(x,y)
    integer,intent(out) :: x,y
    x = 3
    y = 4
  end subroutine get_dim_node_xyz
  
  subroutine get_node_xyz(node_xyz)
    use arrays, only:dp
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
  end subroutine get_node_xyz
  
  subroutine write_elem_file(ELEMFILE)
    use other_consts, only: MAX_FILENAME_LEN
    character(len=MAX_FILENAME_LEN), intent(out) :: ELEMFILE  
    ELEMFILE = "test.ipelem"   
    open(10, file=ELEMFILE, status="replace")
    write(10,*) 'CMISS Version 2.1  ipelem File Version 2'
    write(10,*) 'Heading:'
    write(10,*) ''
    write(10,*) 'The number of elements is [1]: 3'
    write(10,*) ''
    write(10,*) 'Element number [    1]:  1'
    write(10,*) 'The number of geometric Xj-coordinates is [3]: 3'
    write(10,*) 'The basis function type for geometric variable 1 is [1]:  1'
    write(10,*) 'The basis function type for geometric variable 2 is [1]:  1'
    write(10,*) 'The basis function type for geometric variable 3 is [1]:  1'
    write(10,*) 'Enter the 2 global numbers for basis 1:  1  2'
    write(10,*) ''
    write(10,*) 'Element number [ 1002]:  2'
    write(10,*) 'The number of geometric Xj-coordinates is [3]: 3'
    write(10,*) 'The basis function type for geometric variable 1 is [1]:  1'
    write(10,*) 'The basis function type for geometric variable 2 is [1]:  1'
    write(10,*) 'The basis function type for geometric variable 3 is [1]:  1'
    write(10,*) 'Enter the 2 global numbers for basis 1:  2  3'
    write(10,*) ''
    write(10,*) 'Element number [ 1003]:  3'
    write(10,*) 'The number of geometric Xj-coordinates is [3]: 3'
    write(10,*) 'The basis function type for geometric variable 1 is [1]:  1'
    write(10,*) 'The basis function type for geometric variable 2 is [1]:  1'
    write(10,*) 'The basis function type for geometric variable 3 is [1]:  1'
    write(10,*) 'Enter the 2 global numbers for basis 1:  2  4'
    close(10)  
  end subroutine write_elem_file
  
  subroutine delete_elem_file(ELEMFILE)
    use other_consts, only: MAX_FILENAME_LEN
    character(len=MAX_FILENAME_LEN), intent(out) :: ELEMFILE  
    ELEMFILE = "test.ipelem"   
    open(10,file=ELEMFILE, status="old")
    close(10, status="delete") 
  end subroutine delete_elem_file
  
  subroutine get_num_elems(num_elems)
    integer,intent(out) :: num_elems
    num_elems = 3
  end subroutine get_num_elems
  
  subroutine get_dim_elem_nodes(x)
    integer,intent(out) :: x
    call get_num_elems(x)  
  end subroutine get_dim_elem_nodes

  subroutine get_elem_nodes(elem_nodes)
     integer,intent(out) :: elem_nodes(2,3)  
     elem_nodes(1,1)=1
     elem_nodes(2,1)=2
     elem_nodes(1,2)=2
     elem_nodes(2,2)=3
     elem_nodes(1,3)=2
     elem_nodes(2,3)=4
  end subroutine get_elem_nodes
  
  subroutine get_dim_elems_at_node(tree_type,x)
     character(len=20),intent(in) :: tree_type
     integer,intent(out) :: x
     call get_num_nodes(tree_type,x)
  end subroutine get_dim_elems_at_node
  
  subroutine get_elems_at_node(elems_at_node)
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
  end subroutine get_elems_at_node
  
  subroutine get_dim_elem_cnct(x)
    integer, intent(out) :: x
    call get_num_elems(x)  
  end subroutine get_dim_elem_cnct
    
  subroutine get_elem_cnct(elem_cnct)
    integer,intent(out) :: elem_cnct(-1:1,0:2,0:3) !array size: -1:1,0:2,0:num_elems
    elem_cnct = 0
    elem_cnct(1,0,1)=2
    elem_cnct(1,1,1)=2
    elem_cnct(1,2,1)=3
    elem_cnct(-1,0,2)=1
    elem_cnct(-1,1,2)=1
    elem_cnct(-1,0,3)=1
    elem_cnct(-1,1,3)=1      
  end subroutine get_elem_cnct    
  
  subroutine get_dim_elem_field(x)
    integer, intent(out) :: x
    call get_num_elems(x)    
  end subroutine get_dim_elem_field
  
  subroutine get_elem_field(elem_field)
    use indices, only: ne_length,num_ne
    use arrays, only:dp
    real(dp),intent(out) :: elem_field(num_ne,3) ! array size: num_ne, num_elems
    elem_field = 0.0_dp
    elem_field(ne_length,1) = 100.0_dp
    elem_field(ne_length,2) = 70.710678118654755_dp
    elem_field(ne_length,3) = 70.710678118654755_dp  
  end subroutine get_elem_field
  
  subroutine get_dim_elem_direction(x)
    integer, intent(out) :: x
    call get_num_elems(x)     
  end subroutine get_dim_elem_direction
   
  subroutine get_elem_direction(elem_direction)
    use arrays, only:dp
    real(dp),intent(out) :: elem_direction(3,3) !array size: 3,num_elems  
    elem_direction(1,1) = 0.0_dp
    elem_direction(2,1) = 0.0_dp
    elem_direction(3,1) = -1.0_dp
    elem_direction(1,2) = 0.0_dp
    elem_direction(2,2 )= -0.70710678118654746_dp
    elem_direction(3,2) = -0.70710678118654746_dp
    elem_direction(1,3) = 0.0_dp
    elem_direction(2,3) = 0.70710678118654746_dp
    elem_direction(3,3) = -0.70710678118654746_dp  
  end subroutine get_elem_direction

  subroutine get_dim_elem_ordrs(x)
    integer, intent(out) :: x
    call get_num_elems(x)     
  end subroutine get_dim_elem_ordrs

  subroutine get_elem_ordrs(elem_ordrs)
    use indices, only: num_ord
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
  end subroutine get_elem_ordrs
  
  subroutine get_maxgen(maxgen)
    integer,intent(out) :: maxgen
    maxgen = 2
  end subroutine get_maxgen
  
  subroutine get_num_units(num_units)
    integer,intent(out) :: num_units
    num_units = 2
  end subroutine get_num_units
   
  subroutine get_dim_units(x)
    integer,intent(out) :: x
    call get_num_units(x)  
  end subroutine get_dim_units
  
  subroutine get_units(units)
    integer,intent(out) :: units(2)
    units(1)=2
    units(2)=3
  end subroutine get_units
       
  subroutine get_dim_elem_units_below(x)
    integer,intent(out) :: x
    call get_num_elems(x)
  end subroutine get_dim_elem_units_below
  
  subroutine get_elem_units_below(elem_units_below)
    integer,intent(out) :: elem_units_below(1:3) ! 1:num_elems
    elem_units_below(1) = 2 
    elem_units_below(2) = 1 
    elem_units_below(3) = 1
  end subroutine get_elem_units_below
  
end module test_data
