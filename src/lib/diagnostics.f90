!> \file
!> \author Merryn Tawhai
!> \brief This module handles diagnostics.
!>
!> \section LICENSE
!>
!>
!> Contributor(s):
!>
!>\Description
!> This module handles diagnostics
module diagnostics

  implicit none
  integer :: diagnostics_level ! level 0 - no diagnostics level 1 - only prints subroutine names, level 2 - prints sub names and contents of variables

  private
  public enter_exit, get_diagnostics_level, set_diagnostics_level

contains

!!!######################################################################

  subroutine enter_exit(sub_name, state)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ENTER_EXIT" :: ENTER_EXIT
    use other_consts, only: MAX_SUBNAME_LEN
    implicit none

    integer,intent(in) :: state
    character(len=MAX_SUBNAME_LEN), intent(in) :: sub_name

	!print *,"diagnostics_level=",diagnostics_level

    if(diagnostics_level.GE.1)then
      if(state.eq.1)then
        write(*,'('' Entering subroutine '',60A,'':'')') sub_name(1:MAX_SUBNAME_LEN)
      else
        write(*,'('' Exiting subroutine '',60A,'':'')') sub_name(1:MAX_SUBNAME_LEN)
      endif
    endif

  end subroutine enter_exit

!!!######################################################################

  subroutine get_diagnostics_level(level)
    implicit none

    integer :: level

    level = diagnostics_level

  end subroutine get_diagnostics_level

!!!######################################################################

  subroutine set_diagnostics_level(level)
  !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"SO_SET_DIAGNOSTICS_LEVEL":: SET_DIAGNOSTICS_LEVEL
    implicit none

    integer, intent(in) :: level

    diagnostics_level = level


  end subroutine set_diagnostics_level


end module diagnostics
