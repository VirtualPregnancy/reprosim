module nutrient_transport_c
  implicit none
!Interfaces
private

contains

!
!> Perfusion indices
  subroutine evaluate_transport_c() bind(C, name="evaluate_transport_c")

    use nutrient_transport, only: evaluate_transport
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_evaluate_transport()
#else
    call evaluate_transport()
#endif

  end subroutine evaluate_transport_c


end module nutrient_transport_c
