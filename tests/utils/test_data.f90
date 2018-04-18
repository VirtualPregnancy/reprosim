module test_data
!*Description:* This module creates the data required by unit tests. 
! It creates files and populates that are expected by each tested subroutine
!
!*Contributor(s):* Monika Byrne

  implicit none

  private
  public write_node_file, delete_node_file, set_nodes, set_node_xyz

contains
  subroutine write_node_file(NODEFILE)  
    use other_consts, only: MAX_FILENAME_LEN
 
    character(len=MAX_FILENAME_LEN), intent(in) :: NODEFILE
    
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
 
    character(len=MAX_FILENAME_LEN), intent(in) :: NODEFILE
  
    open(10,file=NODEFILE, status="old")
    close(10, status="delete")     
  end subroutine delete_node_file
  
  subroutine set_nodes(nodes)
  
    integer, intent(inout) :: nodes(4)
    
    nodes(1)=1
    nodes(2)=2
    nodes(3)=3
    nodes(4)=4  
  
  end subroutine set_nodes
  
  
  subroutine set_node_xyz(node_xyz)
   use arrays, only:dp
   real(dp),intent(inout) :: node_xyz(3,4)
  
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

end module test_data
