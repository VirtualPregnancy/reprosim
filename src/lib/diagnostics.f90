module diagnostics
!*Description:* This module handles diagnostics. There are 3 diagnostic levels available: level 0 - no diagnostics; level 1- prints subroutine names; and level 2- prints subroutine names as well as variables that each subroutine calculates.
!

  implicit none
  integer :: diagnostics_level ! level 0 - no diagnostics level 1 - only prints subroutine names, level 2 - prints sub names and contents of variables

  private
  public enter_exit, get_diagnostics_level, set_diagnostics_level

contains

!!!######################################################################

  subroutine enter_exit(sub_name, state)
  !*Description:* Prints the subroutine name
    use other_consts, only: MAX_SUBNAME_LEN
    implicit none
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ENTER_EXIT" :: ENTER_EXIT
 
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
    !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"SO_GET_DIAGNOSTICS_LEVEL":: GET_DIAGNOSTICS_LEVEL

    integer :: level

    level = diagnostics_level

  end subroutine get_diagnostics_level

!!!######################################################################

  subroutine set_diagnostics_level(level)
  
    implicit none
    !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"SO_SET_DIAGNOSTICS_LEVEL":: SET_DIAGNOSTICS_LEVEL

    integer, intent(in) :: level

    diagnostics_level = level


  end subroutine set_diagnostics_level


end module diagnostics
