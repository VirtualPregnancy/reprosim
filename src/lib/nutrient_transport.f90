module nutrient_transport
!*Description:* This module contains tools that are used to solve systems of equations nutrient transport .
!
! Descriptions for subroutines that are not included in the subroutine:

  implicit none
  !Module parameters

  !Module types

  !Module depvar

  !Interfaces
  private
  public evaluate_transport

contains

subroutine evaluate_transport()
    use diagnostics, only: enter_exit
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_TRANSPORT" :: EVALUATE_TRANSPORT

    character(len=60) :: sub_name

    sub_name = 'evaluate_transport'

    call enter_exit(sub_name,1)

    call enter_exit(sub_name,2)
end subroutine evaluate_transport

end module nutrient_transport
