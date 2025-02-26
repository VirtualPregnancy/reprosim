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
       ne_radius_out,ne_group,ne_Qdot,ne_viscfact,ne_hb,&
       ne_sa,ne_artvol,ne_artsa,ne_veinvol,ne_veinsa, &
       ne_comp
  !Indices for fetal fields
  integer :: num_nj_fetal, njf_press,njf_vol,njf_comp,njf_type,njf_netQ
  integer :: nef_K, nef_L,nef_dQdt
  ! indices for unit_field
  integer :: num_nu,nu_perf,nu_blood_press

  !model type
  character(len=60) :: model_type

public num_ord,no_gen,no_hord,no_sord

public num_nj,nj_aw_press,nj_bv_press

public num_ne,ne_radius,ne_length,ne_vol,&
      ne_resist,ne_radius_in,ne_radius_out,&
      ne_group,ne_Qdot,ne_viscfact,ne_hb,&
      ne_sa,ne_artvol,ne_artsa,ne_veinvol,ne_veinsa,&
      nef_K, nef_L,nef_dQdt,ne_comp


public num_nu,nu_perf,nu_blood_press

public num_nj_fetal, njf_press,njf_vol,njf_comp,njf_type,njf_netQ

public model_type

!dec$ attributes dllexport :: ne_Qdot, ne_radius, ne_length, num_ord
!dec$ attributes dllexport :: num_nu

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
    num_ne=17
    ne_radius=1 !strained average radius over whole element
    ne_radius_in=2 !strained radius into an element
    ne_radius_out=3 !strained radius out of an element
    ne_length=4!length of an elevent
    ne_vol=5!element volume
    ne_Qdot=7 !flow in an element
    ne_resist=8 !resistance of a blood vessel
    ne_group=9!Groups vessels into arteries (field=0), capillaries (field=1),veins(field=2),anastomoses(field=3)
    ne_viscfact = 10
    ne_hb = 11
    ne_sa = 12
    ne_artvol=13 !Waste of memory
    ne_artsa=14 !Waste of memory
    ne_veinvol = 15 !waste of memory
    ne_veinsa = 16  !waste of memory
    ne_comp = 17 !compliance

    nef_K = 13
    nef_L = 14
    nef_dQdt = 15
    !indices for units
    num_nu=2
    nu_perf=1
    nu_blood_press=2

    num_nj_fetal = 5
    njf_press=1 !pressure in compartment
    njf_vol = 2 !volume of compartment
    njf_comp = 3 !compliance associated with compartment
    njf_type =  4 !type of compartment
    njf_netQ=5


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
