module test_data
!*Description:* This module creates the data required by multiple unit tests. 
!
!*Contributor(s):* Monika Byrne

  implicit none

  private
  public write_node_file, delete_node_file, write_elem_file, delete_elem_file
           
contains
  subroutine write_node_file(NODEFILE)  
    !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"SO_WRITE_NODE_FILE" :: WRITE_NODE_FILE
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
    !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"SO_DELETE_NODE_FILE" :: DELETE_NODE_FILE
    use other_consts, only: MAX_FILENAME_LEN 
    character(len=MAX_FILENAME_LEN), intent(out) :: NODEFILE  
    NODEFILE = "test.ipnode"
    open(10,file=NODEFILE, status="old")
    close(10, status="delete")     
  end subroutine delete_node_file
   
  subroutine write_elem_file(ELEMFILE)
    !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"SO_WRITE_ELEM_FILE" :: WRITE_ELEM_FILE
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
    !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"SO_DELETE_ELEM_FILE" :: DELETE_ELEM_FILE
    use other_consts, only: MAX_FILENAME_LEN
    character(len=MAX_FILENAME_LEN), intent(out) :: ELEMFILE  
    ELEMFILE = "test.ipelem"   
    open(10,file=ELEMFILE, status="old")
    close(10, status="delete") 
  end subroutine delete_elem_file
   
end module test_data
