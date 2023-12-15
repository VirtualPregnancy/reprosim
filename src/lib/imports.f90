module imports
!*Brief Description:* This module contains all the subroutines required to
!import fields, previous model results, etc.
!*LICENSE:*
!
!
!
!*Full Description:*
!
  !
  use arrays
  use diagnostics
  use geometry, only: get_final_real
  use indices
  use other_consts

  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private
  public import_exelemfield

contains
!
!##############################################################################
!
!>*import_exelemfield:* This subroutine reads in the content of an exelem field file (1 field)
 subroutine import_exelemfield(FILENAME,field_no)

   character(len=MAX_FILENAME_LEN),intent(in) :: FILENAME
   integer, intent(in) :: field_no
   !local variables
   integer :: ierror,ne,nunit
   character(LEN=132) :: ctemp1,exfile
   real(dp) :: fieldval

   character(len=60) :: sub_name

   sub_name = 'import_exelemfield'
   call enter_exit(sub_name,1)

   open(10, file=FILENAME, status='old')
   ne = 0
   read_elem_field : do !define a do loop name
     !.......read element flow
     read(unit=10, fmt="(a)", iostat=ierror) ctemp1
     if(index(ctemp1, "Values:")> 0) then
       ne = ne+1
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       call get_final_real(ctemp1,fieldval)
       if(fieldval.lt.0.0_dp) fieldval = zero_tol
         elem_field_fetal(field_no,ne) =fieldval! read it in
       end if
       if(ne.ge.num_elems_fetal) exit read_elem_field
     end do read_elem_field

   close(10)

    call enter_exit(sub_name,2)
 end subroutine import_exelemfield

end module imports