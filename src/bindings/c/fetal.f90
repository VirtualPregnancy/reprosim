module fetal_c
  implicit none
!Interfaces
private

contains

!
!> Perfusion fetal
  subroutine fetal_model_c(OUTDIR,filename_len,dt,num_heart_beats,T_beat,T_vs,T_as,T_v_delay,U0RV,EsysRV,EdiaRV,RvRv,&
          U0LV,EsysLV,EdiaLV,RvLV,U0A,V0V,V0A) bind(C, name="fetal_model_c")
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use arrays, only: dp
    use other_consts, only: MAX_FILENAME_LEN
    use fetal, only: fetal_model
    implicit none

    type(c_ptr), value, intent(in) :: OUTDIR
    real(dp), intent(in) :: dt
    integer, intent(in) :: num_heart_beats
    real(dp), intent(in) :: T_beat
    real(dp), intent(in) :: T_vs
    real(dp),  intent(in) :: T_as
    real(dp),  intent(in) :: T_v_delay
    real(dp),  intent(in) :: U0RV
    real(dp),  intent(in) :: EsysRV
    real(dp),  intent(in) :: EdiaRV
    real(dp),  intent(in) :: RvRV
    real(dp),  intent(in) :: U0LV
    real(dp),  intent(in) :: EsysLV
    real(dp),  intent(in) :: EdiaLV
    real(dp),  intent(in) :: RvLV
    real(dp),  intent(in) :: U0A
    real(dp), intent(in) :: V0V
    real(dp), intent(in) :: V0A
    integer,intent(in) :: filename_len
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, OUTDIR, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_fetal_model(filename_f,dt,num_heart_beats,T_beat,T_vs,T_as,T_v_delay,U0RV,EsysRV,EdiaRV,&
            RvRv,U0LV,EsysLV,EdiaLV,RvLV,U0A,V0V,V0A)
#else
    call fetal_model(filename_f,dt,num_heart_beats,T_beat,T_vs,T_as,T_v_delay,U0RV,EsysRV,EdiaRV,RvRv,&
            U0LV,EsysLV,EdiaLV,RvLV,U0A,V0V,V0A)
#endif

  end subroutine fetal_model_c

  !
!> assign_fetal_arrays
  subroutine assign_fetal_arrays_c() bind(C, name="assign_fetal_arrays_c")

    use fetal, only: assign_fetal_arrays
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_assign_fetal_arrays
#else
    call assign_fetal_arrays
#endif

  end subroutine assign_fetal_arrays_c

end module fetal_c