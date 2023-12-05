
module test_diagnostics
  use testdrive, only : new_unittest, unittest_type, error_type, check
  implicit none
  private

  public :: collect_diagnostics

contains

!> Collect all exported unit tests
subroutine collect_diagnostics(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("test_set_and_get", test_set_and_get) &
    ]

end subroutine collect_diagnostics

subroutine test_set_and_get(error)
  use diagnostics, only: get_diagnostics_level, set_diagnostics_level
  implicit none

  type(error_type), allocatable, intent(out) :: error

  integer :: level

  call get_diagnostics_level(level)
  call check(error, 0, level)
  if (allocated(error)) return

  call set_diagnostics_level(1)
  call get_diagnostics_level(level)
  call check(error, 1, level)
  if (allocated(error)) return
  
end subroutine test_set_and_get

end module test_diagnostics
