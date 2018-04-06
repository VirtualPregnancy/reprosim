
module indices_c
  implicit none
!Interfaces
private

contains
  
!
!> Perfusion indices
  subroutine perfusion_indices_c() bind(C, name="perfusion_indices_c")

    use indices, only: perfusion_indices
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_perfusion_indices()
#else
    call perfusion_indices()
#endif

  end subroutine perfusion_indices_c

  function get_ne_radius_c() result(res) bind(C, name="get_ne_radius_c")

    use indices, only: get_ne_radius
    implicit none
    integer :: res

    res = get_ne_radius()

  end function get_ne_radius_c

end module indices_c
