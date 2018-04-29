module other_consts
!
!*Description:* This module contains definition of constants
!
  use arrays, only: dp
  implicit none

  integer, parameter :: MAX_FILENAME_LEN = 255, MAX_STRING_LEN = 100, MAX_SUBNAME_LEN = 60

  real(dp), parameter :: PI = 3.14159265358979_dp

  private
  public MAX_SUBNAME_LEN, MAX_STRING_LEN, MAX_FILENAME_LEN, PI
end module other_consts
