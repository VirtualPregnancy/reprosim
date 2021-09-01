module indices
! 
!*Description:* This module creates indices used by multiple modules in this package.
!
  implicit none
!parameters
  ! indices for elem_ordrs
  integer :: num_ord=3,no_gen=1,no_hord=2,no_sord=3
  ! indices for node_fields
  integer :: num_nj,nj_aw_press,nj_bv_press
  ! indices for elem_field
  integer ::num_ne,ne_radius,ne_length,ne_vol,&
       ne_resist,ne_radius_in,&
       ne_radius_out,ne_group,ne_Qdot
  ! indices for unit_field
  integer :: num_nu,nu_perf,nu_blood_press

  !model type
  character(len=60) :: model_type

public num_ord,no_gen,no_hord,no_sord

public num_nj,nj_aw_press,nj_bv_press

public num_ne,ne_radius,ne_length,ne_vol,&
      ne_resist,ne_radius_in,ne_radius_out,&
      ne_group,ne_Qdot

public num_nu,nu_perf,nu_blood_press

public model_type

!Interfaces
private
public perfusion_indices, get_ne_radius

contains

!######################################################################
!
!> Perfusion indices
  subroutine perfusion_indices

    use diagnostics, only: enter_exit
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_PERFUSION_INDICES" :: PERFUSION_INDICES
    
    character(len=60) :: sub_name

    sub_name = 'perfusion_indices'
    call enter_exit(sub_name,1)

    ! indices for node_field
    num_nj=1
    nj_bv_press=1 !pressure in blood vessel
    ! indices for elem_field
    num_ne=9
    ne_radius=1 !strained average radius over whole element
    ne_radius_in=2 !strained radius into an element
    ne_radius_out=3 !strained radius out of an element
    ne_length=4!length of an elevent
    ne_vol=5!element volume
    ne_Qdot=7 !flow in an element
    ne_resist=8 !resistance of a blood vessel
    ne_group=9!Groups vessels into arteries (field=0), capillaries (field=1),veins(field=2),anastomoses(field=3)
    !indices for units
    num_nu=2
    nu_perf=1
    nu_blood_press=2

     call enter_exit(sub_name,2)
  end subroutine perfusion_indices

  function get_ne_radius() result(res)
  
    use diagnostics, only: enter_exit
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GET_NE_RADIUS" :: GET_NE_RADIUS

    character(len=60) :: sub_name
    integer :: res

    sub_name = 'get_ne_radius'
    call enter_exit(sub_name,1)

    res=ne_radius

    call enter_exit(sub_name,2)
  end function get_ne_radius

end module indices
