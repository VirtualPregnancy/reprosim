module fetal
!*Description:* This module contains fetal models
! Descriptions for subroutines that are not included in the subroutine:

  use arrays
  use diagnostics
  use indices
  use other_consts
  use pressure_resistance_flow, only: calculate_resistance,capillary_resistance
  use imports, only:  import_exelemfield

  implicit none
  !Module parameters

  !Module types

  !Module depvar

  !Interfaces
  private
  public fetal_model
  public assign_fetal_arrays

  !real(dp),parameter,private :: T_beat = 0.43_dp         ! heart beat period (s)
  !real(dp),parameter,private :: T_vs  = 0.215_dp ! Time period of ventricular contraction (s)
  !real(dp), parameter,private :: T_as = 0.1075_dp !Time period of atrial contraction (s)
  !real(dp), parameter,private :: T_v_delay = 0.1075_dp !delay in ventrial contraction (compare to atria) (s)
  !real(dp), parameter, private :: U0RV = 5332.89_dp !Pa
  !real(dp), parameter, private :: EsysRV = 0.399967_dp !Pa/mm3
  !real(dp), parameter, private :: EdiaRV = 0.0399967_dp !Pa/mm3
  !real(dp), parameter, private :: RvRV = 0.002!0.010665 !Pa.s/mm3
  !real(dp), parameter, private :: U0LV = 5332.89_dp !Pa
  !real(dp), parameter, private :: EsysLV = 0.399967_dp !Pa/mm3
  !real(dp), parameter, private :: EdiaLV = 0.0399967_dp !Pa/mm3
  !real(dp), parameter, private :: RvLV = 0.002!0.010665_dp !Pa.s/mm3
  !real(dp), parameter, private :: U0A = 399.967_dp !Pa


contains
    subroutine fetal_model(OUTDIR,dt,num_heart_beats,T_beat,T_vs,T_as,T_v_delay,U0RV,EsysRV,EdiaRV,RvRv,U0LV,EsysLV,EdiaLV,&
            RvLV,U0A,V0V,V0A)
        use diagnostics, only: enter_exit,get_diagnostics_level
        use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_FETAL_MODEL" :: FETAL_MODEL

        character(len=MAX_FILENAME_LEN), intent(in) :: OUTDIR
        real(dp), intent(in) :: dt                     !Time step (s)
        integer, intent(in) :: num_heart_beats         !num heart beats
        real(dp), intent(in) :: T_beat! = 0.43_dp         ! heart beat period (s)
        real(dp), intent(in) :: T_vs ! = 0.215_dp ! Time period of ventricular contraction (s)
        real(dp),  intent(in) :: T_as != 0.1075_dp !Time period of atrial contraction (s)
        real(dp),  intent(in) :: T_v_delay != 0.1075_dp !delay in ventrial contraction (compare to atria) (s)
        real(dp),  intent(in) :: U0RV! = 5332.89_dp !Pa
        real(dp),  intent(in) :: EsysRV != 0.399967_dp !Pa/mm3
        real(dp),  intent(in) :: EdiaRV! = 0.0399967_dp !Pa/mm3
        real(dp),  intent(in) :: RvRV! = 0.002!0.010665 !Pa.s/mm3
        real(dp),  intent(in) :: U0LV! = 5332.89_dp !Pa
        real(dp),  intent(in) :: EsysLV! = 0.399967_dp !Pa/mm3
        real(dp),  intent(in) :: EdiaLV! = 0.0399967_dp !Pa/mm3
        real(dp),  intent(in) :: RvLV! = 0.002!0.010665_dp !Pa.s/mm3
        real(dp),  intent(in) :: U0A! = 399.967_dp !Pa
        real(dp), intent(in) :: V0V ! Initial ventricle volume
        real(dp), intent(in) :: V0A !Initial atrial volume

        real(dp) :: time                  !current time (s)
        real(dp) :: T_interval            ! the total length of the heat beat (s)


        integer :: n    !current heart beat
        real(dp) :: ttime !time within current heart beat
        integer :: np,ne,np_in,np_out

        real(dp) :: Avent !Ventricular activation (no units)
        real(dp) :: Aatria !Atrial activation (no units)
        real(dp) :: dpress,Pgrad,Qnod,dQ,Vnod,press
        real(dp) :: art_resistance,ven_resistance,total_volume
        character(len=60) :: mesh_type
        logical :: continue
        character(len=60) :: sub_name
        character(len=MAX_FILENAME_LEN) :: filename
        integer :: diagnostics_level

        !------
        sub_name = 'assign_fetal_arrays'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        Avent=0.0_dp !initialise
        Aatria = 0.0_dp !Initialise

        T_interval = num_heart_beats * T_beat
        write(*,*) "Simulating for" , T_interval, " s"

        node_field_fetal(njf_type,:) = node_xyz_fetal(1,:)
        node_field_fetal(njf_press,:) = node_xyz_fetal(2,:)
        node_field_fetal(njf_comp,:) = node_xyz_fetal(3,:)
        filename = TRIM(OUTDIR) //  'group.exelem'
        call import_exelemfield(filename,ne_group)
        filename = TRIM(OUTDIR) // 'R.exelem'
        call import_exelemfield(filename,ne_resist)
        filename = TRIM(OUTDIR) // 'K.exelem'
        call import_exelemfield(filename,nef_K)
        filename = TRIM(OUTDIR) // 'L.exelem'
        call import_exelemfield(filename,nef_L)
        total_volume = 0.0_dp
        do np = 1,num_nodes_fetal
            if(node_field_fetal(njf_type,np).le.2.)then
                node_field_fetal(njf_vol,np) = V0V!node_field_fetal(njf_press,np)*node_field_fetal(njf_comp,np)
            elseif(node_field_fetal(njf_type,np).eq.3.)then
                node_field_fetal(njf_vol,np) = V0A!3000.
            else
                node_field_fetal(njf_vol,np) = node_field_fetal(njf_press,np)*node_field_fetal(njf_comp,np)
            endif
            total_volume = total_volume + node_field_fetal(njf_vol,np)
        end do

        write(*,*) 'Total  blood volume (ml):',total_volume/1000.

        write(*,*) 'Calculating placental resistance'
        mesh_type = 'simple_tree'
        elem_field(ne_viscfact,:) = 1.0_dp !initialise viscosity factor
        call calculate_resistance(0.33600e-02_dp,mesh_type)
        call tree_resistance(art_resistance,ven_resistance)
        write(*,*) 'Arterial resistance (Pa.s/mm3)= ', art_resistance
        write(*,*) 'Venous resistance (Pa.s/mm3)= ', ven_resistance
        do ne =1,num_elems_fetal
            if (abs(elem_field_fetal(ne_group,ne)-9.0_dp).lt.loose_tol)then!Umbilical artery
                elem_field_fetal(ne_resist,ne) = art_resistance ! Pa s /mm3
            elseif(abs(elem_field_fetal(ne_group,ne)-10.0_dp).lt.loose_tol)then!Umbilical vein
                elem_field_fetal(ne_resist,ne) = ven_resistance ! Pa s /mm3
            end if
        end do

        Write(*,*) 'Initialising flows'
        !Initialise flows
        do ne =1,num_elems_fetal
           np_in = elem_nodes_fetal(1,ne)
           np_out = elem_nodes_fetal(2,ne)
           Pgrad = node_field_fetal(njf_press,np_in)-node_field_fetal(njf_press,np_out)
            if (abs(elem_field_fetal(ne_group,ne)-2.0_dp).lt.loose_tol)then!R-Q unit
               call rq_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne))
            else if (abs(elem_field_fetal(ne_group,ne)-9.0_dp).lt.loose_tol)then!R-Q unit
               call rq_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne))
            else if (abs(elem_field_fetal(ne_group,ne)-10.0_dp).lt.loose_tol)then!R-Q unit
               call rq_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne))
            else if (abs(elem_field_fetal(ne_group,ne)-3.0_dp).lt.loose_tol)then!R-Q-L units, initialise to zero L and will be time stepped units
                call rq_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne))
            elseif (abs(elem_field_fetal(ne_group,ne)-4.0_dp).lt.loose_tol)then!FOROMAN OVALE
               if(Pgrad.le.0.0_dp)then
                   elem_field_fetal(ne_Qdot,ne)=0.0_dp
               else
                   elem_field_fetal(ne_Qdot,ne)=(Pgrad/elem_field_fetal(nef_k,ne))**(1.0_dp/0.625_dp)
               end if
             elseif  (abs(elem_field_fetal(ne_group,ne)-5.0_dp).lt.loose_tol)then!R-Q-K unit (ductus venosus
                call rqk_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne),&
                           elem_field_fetal(nef_K,ne))
             elseif  (abs(elem_field_fetal(ne_group,ne)-6.0_dp).lt.loose_tol)then!R-Q-K-L unit (ductus arterious, will time step
                call rqk_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne),&
                           elem_field_fetal(nef_K,ne))
             elseif  (abs(elem_field_fetal(ne_group,ne)-8.0_dp).lt.loose_tol)then!Cardiac exit
                if(Pgrad.lt.0.0_dp)then
                    elem_field_fetal(ne_Qdot,ne)=0.0_dp
                else
                    call rqk_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne),&
                           elem_field_fetal(nef_K,ne))
                end if
            else
            end if
        end do

        time = 0.0_dp !initialise the simulation time.
        ttime = 0.0_dp !initialise the simulation time.
        open(10, file=TRIM(OUTDIR) // 'results_volume.out', status='replace')
        open(20, file=TRIM(OUTDIR) // 'results_pressure.out', status='replace')
        open(30, file=TRIM(OUTDIR) // 'results_flow.out', status='replace')
        open(40, file=TRIM(OUTDIR) // 'results_element_flow.out', status='replace')
        WRITE(10,'(26(F15.4,X,","),(F15.4,X))')&
                time, ttime,Avent,Aatria,node_field_fetal(njf_vol,1),&
                        node_field_fetal(njf_vol,2), node_field_fetal(njf_vol,3),node_field_fetal(njf_vol,4),&
                node_field_fetal(njf_vol,5), node_field_fetal(njf_vol,6),node_field_fetal(njf_vol,7),&
                node_field_fetal(njf_vol,8), node_field_fetal(njf_vol,9),node_field_fetal(njf_vol,10),&
                node_field_fetal(njf_vol,11), node_field_fetal(njf_vol,12),node_field_fetal(njf_vol,13),&
                node_field_fetal(njf_vol,14), node_field_fetal(njf_vol,15),node_field_fetal(njf_vol,16),&
                node_field_fetal(njf_vol,17), node_field_fetal(njf_vol,18),node_field_fetal(njf_vol,19),&
                node_field_fetal(njf_vol,20), node_field_fetal(njf_vol,21),node_field_fetal(njf_vol,22),&
                node_field_fetal(njf_vol,23)

        WRITE(20,'(26(F15.4,X,","),(F15.4,X))')&
                time, ttime,Avent,Aatria,node_field_fetal(njf_press,1),&
                        node_field_fetal(njf_press,2), node_field_fetal(njf_press,3),node_field_fetal(njf_press,4),&
                node_field_fetal(njf_press,5), node_field_fetal(njf_press,6),node_field_fetal(njf_press,7),&
                node_field_fetal(njf_press,8), node_field_fetal(njf_press,9),node_field_fetal(njf_press,10),&
                node_field_fetal(njf_press,11), node_field_fetal(njf_press,12),node_field_fetal(njf_press,13),&
                node_field_fetal(njf_press,14), node_field_fetal(njf_press,15),node_field_fetal(njf_press,16),&
                node_field_fetal(njf_press,17), node_field_fetal(njf_press,18),node_field_fetal(njf_press,19),&
                node_field_fetal(njf_press,20), node_field_fetal(njf_press,21),node_field_fetal(njf_press,22),&
                node_field_fetal(njf_press,23)

        WRITE(30,'(26(F15.4,X,","),(F15.4,X))')&
                time, ttime,Avent,Aatria,node_field_fetal(njf_netQ,1),&
                        node_field_fetal(njf_netQ,2), node_field_fetal(njf_netQ,3),node_field_fetal(njf_netQ,4),&
                node_field_fetal(njf_netQ,5), node_field_fetal(njf_netQ,6),node_field_fetal(njf_netQ,7),&
                node_field_fetal(njf_netQ,8), node_field_fetal(njf_netQ,9),node_field_fetal(njf_netQ,10),&
                node_field_fetal(njf_netQ,11), node_field_fetal(njf_netQ,12),node_field_fetal(njf_netQ,13),&
                node_field_fetal(njf_netQ,14), node_field_fetal(njf_netQ,15),node_field_fetal(njf_netQ,16),&
                node_field_fetal(njf_netQ,17), node_field_fetal(njf_netQ,18),node_field_fetal(njf_netQ,19),&
                node_field_fetal(njf_netQ,20), node_field_fetal(njf_netQ,21),node_field_fetal(njf_netQ,22),&
                node_field_fetal(njf_netQ,23)



        continue = .true.
        n = 0
        do while (continue)
            n = n + 1 ! increment the heart beat number
            ttime = 0.0_dp ! each breath starts with ttime=0
            write(*,*) 'Initiating beat number',n
            do while (ttime.lt.T_beat)
                ttime = ttime + dt ! increment the heartbeat time
                time = time + dt ! increment the whole simulation time
                if ((ttime.ge.T_v_delay).and.(ttime.le.T_vs+T_v_delay)) then
                    Avent = sin(pi/T_vs*(ttime-T_v_delay))
                else
                    Avent = 0.0_dp
                end if
                if ((ttime.le.T_as))then
                    Aatria = sin(pi/T_as*(ttime))
                else
                    Aatria = 0.0_dp
                end if
                do np =1,num_nodes_fetal
                    Qnod = 0.0_dp !Net flow passing through the node
                    do ne = 1, elems_at_node_fetal(np,0)
                        if (elem_nodes_fetal(1,elems_at_node_fetal(np,ne)).eq.np)then
                           Qnod = Qnod-elem_field_fetal(ne_Qdot,elems_at_node_fetal(np,ne))
                        else
                            Qnod = Qnod+elem_field_fetal(ne_Qdot,elems_at_node_fetal(np,ne))
                        end if
                    end do !calculating net flow through the nodes
                    node_field_fetal(njf_netQ,np)=Qnod
                    Vnod = node_field_fetal(njf_vol,np)+dt*Qnod
                    node_field_fetal(njf_vol,np) = Vnod
                    if(abs(node_field_fetal(njf_type,np)-1.0_dp).lt.loose_tol)then!right ventricle
                        call ventricle_pressure_step(press,Avent,U0RV,EdiaRV,EsysRV,RvRV,Qnod,Vnod)
                        node_field_fetal(njf_press,np) = press
                    elseif(abs(node_field_fetal(njf_type,np)-2.0_dp).lt.loose_tol)then!left ventricle)
                        call ventricle_pressure_step(press,Avent,U0LV,EdiaLV,EsysLV,RvLV,Qnod,Vnod)
                        node_field_fetal(njf_press,np) = press
                    elseif(abs(node_field_fetal(njf_type,np)-3.0_dp).lt.loose_tol)then!Its an atrium
                        call atrium_pressure_step(press,dt,Aatria,U0A,node_field_fetal(njf_comp,np),Qnod,dQ,Vnod)
                        node_field_fetal(njf_press,np) = press
                    else ! This is a standard node
                        call compartment_pressure_step(dpress,dt,node_field_fetal(njf_comp,np), Qnod,Vnod)
                        node_field_fetal(njf_press,np) = node_field_fetal(njf_press,np) + dpress
                    end if !what type of compartment
                end do !nodal balances

                !! Update flows
                do ne =1,num_elems_fetal !fetal elements
                    np_in = elem_nodes_fetal(1,ne)
                    np_out = elem_nodes_fetal(2,ne)
                    Pgrad = node_field_fetal(njf_press,np_in)-node_field_fetal(njf_press,np_out)
                    if (abs(elem_field_fetal(ne_group,ne)-2.0_dp).lt.loose_tol)then!R-Q unit
                       call rq_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne))
                    else if (abs(elem_field_fetal(ne_group,ne)-9.0_dp).lt.loose_tol)then!R-Q unit
                       call rq_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne))
                    else if (abs(elem_field_fetal(ne_group,ne)-10.0_dp).lt.loose_tol)then!R-Q unit
                       call rq_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne))
                    else if (abs(elem_field_fetal(ne_group,ne)-3.0_dp).lt.loose_tol)then!R-Q-L units, initialise to zero L and will be time stepped units
                        call rql_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne),&
                            elem_field_fetal(nef_L,ne))
                    elseif (abs(elem_field_fetal(ne_group,ne)-4.0_dp).lt.loose_tol)then!FOROMAN OVALE
                        if(Pgrad.le.0.0_dp)then
                           elem_field_fetal(ne_Qdot,ne)=0.0_dp
                        else
                           elem_field_fetal(ne_Qdot,ne)=(Pgrad/elem_field_fetal(nef_k,ne))**(1.0_dp/0.625_dp)
                        end if
                    elseif  (abs(elem_field_fetal(ne_group,ne)-5.0_dp).lt.loose_tol)then!R-Q-K unit (ductus venosus
                        call rqk_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne),&
                           elem_field_fetal(nef_K,ne))
                     elseif  (abs(elem_field_fetal(ne_group,ne)-6.0_dp).lt.loose_tol)then!R-Q-K-L unit (ductus arterious, will time step
                        call rqkl_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne),&
                               elem_field_fetal(nef_K,ne),elem_field_fetal(nef_L,ne))
                     elseif  (abs(elem_field_fetal(ne_group,ne)-8.0_dp).lt.loose_tol)then!Cardiac exits
                        if(Pgrad.lt.0.0_dp)then
                            elem_field_fetal(ne_Qdot,ne)=0.0_dp
                        else
                            call rqk_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne),&
                               elem_field_fetal(nef_K,ne))
                        end if
                    elseif  (abs(elem_field_fetal(ne_group,ne)-1.0_dp).lt.loose_tol)then!Cardiac valves
                        if(Pgrad.lt.0.0_dp)then
                            elem_field_fetal(ne_Qdot,ne)=0.0_dp
                        else
                            call rqkl_unit(dt,elem_field_fetal(ne_Qdot,ne),Pgrad,elem_field_fetal(ne_resist,ne),&
                                elem_field_fetal(nef_K,ne),elem_field_fetal(nef_L,ne))
                        end if
                    end if
                enddo
                WRITE(40,'(33(F15.4,X,","),(F15.4,X))')&
                time, ttime,elem_field_fetal(ne_Qdot,1),elem_field_fetal(ne_Qdot,2),elem_field_fetal(ne_Qdot,3),&
                        elem_field_fetal(ne_Qdot,4),elem_field_fetal(ne_Qdot,5),elem_field_fetal(ne_Qdot,6),&
                        elem_field_fetal(ne_Qdot,7),elem_field_fetal(ne_Qdot,8),elem_field_fetal(ne_Qdot,9),&
                        elem_field_fetal(ne_Qdot,10),elem_field_fetal(ne_Qdot,11),elem_field_fetal(ne_Qdot,12),&
                        elem_field_fetal(ne_Qdot,13),elem_field_fetal(ne_Qdot,14),elem_field_fetal(ne_Qdot,15),&
                        elem_field_fetal(ne_Qdot,16),elem_field_fetal(ne_Qdot,17),elem_field_fetal(ne_Qdot,18),&
                        elem_field_fetal(ne_Qdot,19),elem_field_fetal(ne_Qdot,20),elem_field_fetal(ne_Qdot,21),&
                        elem_field_fetal(ne_Qdot,22),elem_field_fetal(ne_Qdot,23),elem_field_fetal(ne_Qdot,24),&
                        elem_field_fetal(ne_Qdot,25),elem_field_fetal(ne_Qdot,26),elem_field_fetal(ne_Qdot,27),&
                        elem_field_fetal(ne_Qdot,28),elem_field_fetal(ne_Qdot,29),elem_field_fetal(ne_Qdot,20),&
                        elem_field_fetal(ne_Qdot,31),elem_field_fetal(ne_Qdot,32)

                WRITE(10,'(26(F15.4,X,","),(F15.4,X))')&
                time, ttime,Avent,Aatria,node_field_fetal(njf_vol,1),&
                        node_field_fetal(njf_vol,2), node_field_fetal(njf_vol,3),node_field_fetal(njf_vol,4),&
                node_field_fetal(njf_vol,5), node_field_fetal(njf_vol,6),node_field_fetal(njf_vol,7),&
                node_field_fetal(njf_vol,8), node_field_fetal(njf_vol,9),node_field_fetal(njf_vol,10),&
                node_field_fetal(njf_vol,11), node_field_fetal(njf_vol,12),node_field_fetal(njf_vol,13),&
                node_field_fetal(njf_vol,14), node_field_fetal(njf_vol,15),node_field_fetal(njf_vol,16),&
                node_field_fetal(njf_vol,17), node_field_fetal(njf_vol,18),node_field_fetal(njf_vol,19),&
                node_field_fetal(njf_vol,20), node_field_fetal(njf_vol,21),node_field_fetal(njf_vol,22),&
                node_field_fetal(njf_vol,23)
                WRITE(20,'(26(F15.4,X,","),(F15.4,X))')&
                time, ttime,Avent,Aatria,node_field_fetal(njf_press,1),&
                        node_field_fetal(njf_press,2), node_field_fetal(njf_press,3),node_field_fetal(njf_press,4),&
                node_field_fetal(njf_press,5), node_field_fetal(njf_press,6),node_field_fetal(njf_press,7),&
                node_field_fetal(njf_press,8), node_field_fetal(njf_press,9),node_field_fetal(njf_press,10),&
                node_field_fetal(njf_press,11), node_field_fetal(njf_press,12),node_field_fetal(njf_press,13),&
                node_field_fetal(njf_press,14), node_field_fetal(njf_press,15),node_field_fetal(njf_press,16),&
                node_field_fetal(njf_press,17), node_field_fetal(njf_press,18),node_field_fetal(njf_press,19),&
                node_field_fetal(njf_press,20), node_field_fetal(njf_press,21),node_field_fetal(njf_press,22),&
                node_field_fetal(njf_press,23)

                WRITE(30,'(26(F15.4,X,","),(F15.4,X))')&
                time, ttime,Avent,Aatria,node_field_fetal(njf_netQ,1),&
                        node_field_fetal(njf_netQ,2), node_field_fetal(njf_netQ,3),node_field_fetal(njf_netQ,4),&
                node_field_fetal(njf_netQ,5), node_field_fetal(njf_netQ,6),node_field_fetal(njf_netQ,7),&
                node_field_fetal(njf_netQ,8), node_field_fetal(njf_netQ,9),node_field_fetal(njf_netQ,10),&
                node_field_fetal(njf_netQ,11), node_field_fetal(njf_netQ,12),node_field_fetal(njf_netQ,13),&
                node_field_fetal(njf_netQ,14), node_field_fetal(njf_netQ,15),node_field_fetal(njf_netQ,16),&
                node_field_fetal(njf_netQ,17), node_field_fetal(njf_netQ,18),node_field_fetal(njf_netQ,19),&
                node_field_fetal(njf_netQ,20), node_field_fetal(njf_netQ,21),node_field_fetal(njf_netQ,22),&
                node_field_fetal(njf_netQ,23)
            end do
            if (n.eq.num_heart_beats) then
                continue = .false.
            endif
        end do !do while
        close(10)
        close(20)
        close(30)
        close(40)

        total_volume = 0.0_dp
        do np = 1,num_nodes_fetal
            total_volume = total_volume + node_field_fetal(njf_vol,np)
        end do

        write(*,*) 'Total  blood volume (ml):',total_volume/1000.

        call deallocate_fetal_memory

        call enter_exit(sub_name,2)
    end subroutine fetal_model


    subroutine assign_fetal_arrays

        use arrays,only: dp,elem_field_fetal, num_elems_fetal, elem_nodes_fetal, nodes_fetal, elems_fetal, num_nodes_fetal,&
          node_field_fetal,node_xyz_fetal, elem_cnct_fetal,elem_direction_fetal,elems_at_node_fetal,elem_field, num_elems,&
          elem_nodes, nodes, elems, num_nodes, node_field, elem_cnct,elem_direction,elems_at_node,node_xyz
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ASSIGN_FETAL_ARRAYS" :: ASSIGN_FETAL_ARRAYS

        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'assign_fetal_arrays'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        num_nodes_fetal = num_nodes
        if(allocated(nodes_fetal)) deallocate (nodes_fetal)
        allocate (nodes_fetal(num_nodes_fetal))
        nodes_fetal = nodes
        if(allocated(node_xyz_fetal)) deallocate (node_xyz_fetal)
        allocate (node_xyz_fetal(3,num_nodes_fetal))
        node_xyz_fetal = node_xyz
        if(allocated(node_field_fetal)) deallocate (node_field_fetal)
        allocate (node_field_fetal(num_nj_fetal,num_nodes_fetal))
        node_field_fetal = 0.0_dp
        num_elems_fetal = num_elems
        if(allocated(elems_fetal)) deallocate(elems_fetal) !Array that defines nodal connections between elements
        allocate(elems_fetal(num_elems_fetal))
        elems_fetal = elems
        if(allocated(elem_cnct_fetal)) deallocate(elem_cnct_fetal) !Array that defines connections between elements
        allocate(elem_cnct_fetal(-1:1,0:10,0:num_elems_fetal))!Allows up to 10 elements per node
        elem_cnct_fetal = elem_cnct
        if(allocated(elem_nodes_fetal)) deallocate(elem_nodes_fetal)
        allocate(elem_nodes_fetal(2,num_elems_fetal)) !defines in and out nodes at each element
        elem_nodes_fetal = elem_nodes
        if(allocated(elems_at_node_fetal)) deallocate(elems_at_node_fetal)
        allocate(elems_at_node_fetal(num_nodes_fetal,0:10)) !Allows up to 10 elements per node
        elems_at_node_fetal = elems_at_node
        if(allocated(elem_field_fetal)) deallocate(elem_field_fetal)
        allocate(elem_field_fetal(num_ne,num_elems_fetal))
        elem_field_fetal = elem_field
        elem_field_fetal = 0.0_dp
        if(allocated(elem_direction_fetal)) deallocate(elem_direction_fetal)
        allocate(elem_direction_fetal(3,num_elems_fetal))
        elem_direction_fetal = elem_direction



        if(allocated(nodes)) deallocate (nodes)
        if(allocated(node_xyz)) deallocate (node_xyz)
        if(allocated(node_field)) deallocate (node_field)
        if(allocated(elems)) deallocate(elems)
        if(allocated(elem_cnct)) deallocate(elem_cnct)
        if(allocated(elem_nodes)) deallocate(elem_nodes)
        if(allocated(elems_at_node)) deallocate(elems_at_node)
        if(allocated(elem_field)) deallocate(elem_field)
        if(allocated(elem_direction)) deallocate(elem_direction)

        call enter_exit(sub_name,2)

    end subroutine assign_fetal_arrays

    subroutine deallocate_fetal_memory

        use arrays,only: dp,elem_field_fetal, num_elems_fetal, elem_nodes_fetal, nodes_fetal, elems_fetal, num_nodes_fetal,&
          node_field_fetal,node_xyz_fetal, elem_cnct_fetal,elem_direction_fetal,elems_at_node_fetal,elem_field, num_elems,&
          elem_nodes, nodes, elems, num_nodes, node_field, elem_cnct,elem_direction,elems_at_node,node_xyz
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEALLOCATE_FETAL_MEMORY" :: DEALLOCATE_FETAL_MEMORY

        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'deallocate_fetal_memory'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        if(allocated(nodes_fetal)) deallocate (nodes_fetal)
        if(allocated(node_xyz_fetal)) deallocate (node_xyz_fetal)
        if(allocated(node_field_fetal)) deallocate (node_field_fetal)
        if(allocated(elems_fetal)) deallocate(elems_fetal) !Array that defines nodal connections between elements
        if(allocated(elem_cnct_fetal)) deallocate(elem_cnct_fetal) !Array that defines connections between elements
        if(allocated(elem_nodes_fetal)) deallocate(elem_nodes_fetal)
        if(allocated(elems_at_node_fetal)) deallocate(elems_at_node_fetal)
        if(allocated(elem_field_fetal)) deallocate(elem_field_fetal)
        if(allocated(elem_direction_fetal)) deallocate(elem_direction_fetal)


        call enter_exit(sub_name,2)

    end subroutine deallocate_fetal_memory

    subroutine ventricle_pressure_step(press,Avent,U0,Edia,Esys,Rv,Q,V)
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_VENTRICLE_PRESSURE_STEP" :: VENTRICLE_PRESSURE_STEP
        real(dp), intent(out) :: press
        real(dp), intent(in) :: Avent
        real(dp), intent(in) :: U0
        real(dp), intent(in) :: Edia
        real(dp), intent(in) :: Esys
        real(dp), intent(in) :: Rv
        real(dp), intent(in) :: Q
        real(dp), intent(in) :: V

        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'ventricle pressure step'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        press = U0*Avent + (Edia + Esys*Avent)*V + Rv*Q

        call enter_exit(sub_name,2)

    end subroutine ventricle_pressure_step

subroutine atrium_pressure_step(dpress,dt,Aatria,U0,comp,Q,dQ,V)
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ATRIUM_PRESSURE_STEP" :: ATRIUM_PRESSURE_STEP
        real(dp), intent(out) :: dpress
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: Aatria
        real(dp), intent(in) :: U0
        real(dp), intent(in) :: comp
        real(dp), intent(in) :: Q
        real(dp), intent(in) :: dQ
        real(dp), intent(in) :: V
        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'atrium pressure step'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        dpress = (U0*Aatria + V/(comp))
        call enter_exit(sub_name,2)

    end subroutine atrium_pressure_step


subroutine compartment_pressure_step(dpress,dt,comp,Q,V)
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_COMPARTMENT_PRESSURE_STEP" :: COMPARTMENT_PRESSURE_STEP
        real(dp), intent(out) :: dpress
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: comp
        real(dp), intent(in) :: Q
        real(dp), intent(in) :: V

        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'compartment pressure step'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        dpress = dt*Q/comp
        call enter_exit(sub_name,2)

    end subroutine compartment_pressure_step

!subroutine one_way_valve(dt,dQ,Q,Pgrad, R, K, L)
!        use diagnostics, only: enter_exit,get_diagnostics_level
!
!    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ONE_WAY_VALVE :: ONE_WAY_VALVE
!        real(dp), intent(in) :: dt
!        real(dp), intent(inout) :: dQ
!        real(dp), intent(inout) :: Q
!        real(dp), intent(in) :: Pgrad
!        real(dp), intent(in) :: R
!        real(dp), intent(in) :: K
!        real(dp), intent(in) :: L
!
!        real(dp) :: check_sign
!        character(len=60) :: sub_name
!        integer :: diagnostics_level
!
!        !------
!        sub_name = 'one_way_valve'
!        call enter_exit(sub_name,1)
!        call get_diagnostics_level(diagnostics_level)
!        if(L.gt.0)then
!            dQ = dt*(Pgrad-K*Q**2.0_dp)/L
!            Q = Q+dQ
!        else
!            if(Pgrad.gt.0) then
!                Q = sqrt(Pgrad/K)
!            else
!                Q = -1.0_dp* sqrt(abs(Pgrad)/K)
!            end if
!        end if
!        if(Q.lt.0.0_dp)then
!            Q=0.0_dp
!        end if
!        call enter_exit(sub_name,2)
!    end subroutine one_way_valve

    subroutine rq_unit(dt,Q,Pgrad, R)
        use diagnostics, only: enter_exit,get_diagnostics_level

        !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_RQ_UNIT" :: RQ_UNIT
        real(dp), intent(in) :: dt
        real(dp), intent(inout) :: Q
        real(dp), intent(in) :: Pgrad
        real(dp), intent(in) :: R


        real(dp) :: check_sign
        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'rq_unit'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        Q = Pgrad/R

        call enter_exit(sub_name,2)

    end subroutine rq_unit

        subroutine rql_unit(dt,Q,Pgrad, R,L)
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_RQL_unit" :: RQL_UNIT
        real(dp), intent(in) :: dt
        real(dp), intent(inout) :: Q
        real(dp), intent(in) :: Pgrad
        real(dp), intent(in) :: R
        real(dp), intent(in) :: L

        real(dp) :: dQ
        real(dp) :: check_sign
        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'rql_unit'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        dQ = dt*(Pgrad - R*Q)/L

        Q = Q+dQ


        call enter_exit(sub_name,2)

    end subroutine rql_unit
    subroutine rqk_unit(dt,Q,Pgrad, R,K)
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_RQK_UNIT" :: RQK_UNIT
        real(dp), intent(in) :: dt
        real(dp), intent(inout) :: Q
        real(dp), intent(in) :: Pgrad
        real(dp), intent(in) :: R
        real(dp), intent(in) :: K
        !Subroutine for R-Q-K unit beta = 2


        real(dp) :: check_sign
        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'rqk_unit'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        check_sign = R**2.0_dP + 4.0_dp*K*Pgrad
        if(check_sign.le.0.0_dp)then
            Q=0.0_dp
        else
            Q = (-R + sqrt(check_sign))/(2.0_dp*K) !Taking the largest real root
        end if

        call enter_exit(sub_name,2)

        end subroutine rqk_unit

       subroutine rqkl_unit(dt,Q,Pgrad, R,K,L)
        use diagnostics, only: enter_exit,get_diagnostics_level

    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_RQKL_UNIT" :: RQKL_UNIT
        real(dp), intent(in) :: dt
        real(dp), intent(inout) :: Q
        real(dp), intent(in) :: Pgrad
        real(dp), intent(in) :: R
        real(dp), intent(in) :: K
        real(dp), intent(in) :: L


        real(dp) :: dQ
        character(len=60) :: sub_name
        integer :: diagnostics_level

        !------
        sub_name = 'rqkl_unit'
        call enter_exit(sub_name,1)
        call get_diagnostics_level(diagnostics_level)

        dQ = dt*(Pgrad - k*Q**2.0_dp - R*Q)/L
        Q = Q+dQ


        call enter_exit(sub_name,2)

    end subroutine rqkl_unit

subroutine calculate_compliance(E,Ptm)
    use arrays,only: dp,num_nodes,num_elems,elem_nodes,elem_field
    use indices
    use diagnostics, only: enter_exit,get_diagnostics_level

    real(dp), intent(in) :: E
    real(dp), intent(in) :: Ptm
    character(len=60) :: sub_name
    integer :: diagnostics_level

    !local variables
    integer :: ne,nn,np,ny
    real(dp) :: R0,old_vol,new_vol,h,Rnew
    real(dp) :: total_vol_old,total_vol_new,compliance,cap_res
    character(len=60) :: vessel_type,rheol_type

  sub_name = 'calculate_compliance'
  call enter_exit(sub_name,1)
  call get_diagnostics_level(diagnostics_level)

  vessel_type = 'elastic'
  rheol_type = 'constant_visc'
  total_vol_old = 0.0_dp
  total_vol_new = 0.0_dp
  do ne=1,num_elems
        R0 = elem_field(ne_radius,ne)
        if(R0.gt.0.125_dp)then
            h=0.2_dp*R0
        else
            h=0.8_dp*R0
        endif
        Rnew = R0+3.0_dp*R0**2*Ptm/(4.0_dp*E*h)
        elem_field(ne_comp,ne) = (pi*Rnew**2.0_dp*elem_field(ne_length,ne) - &
                pi*R0**2.0_dp*elem_field(ne_length,ne))/Ptm !mm3/Pa

        total_vol_old = total_vol_old + pi*R0**2.0_dp*elem_field(ne_length,ne)
        total_vol_new = total_vol_new + pi*Rnew**2.0_dp*elem_field(ne_length,ne)
        !add umbilical cord
        if(ne.eq.1) then
            total_vol_old = total_vol_old + pi*R0**2.0_dp*50.0_dp
            total_vol_new = total_vol_new + pi*Rnew**2.0_dp*50.0_dp
        else if (ne.eq.num_arterial_elems+1) then
            total_vol_old = total_vol_old + pi*R0**2.0_dp*50.0_dp
            total_vol_new = total_vol_new + pi*Rnew**2.0_dp*50.0_dp
        end if
  enddo
    compliance = (total_vol_new-total_vol_old)/Ptm

  call enter_exit(sub_name,2)
end subroutine calculate_compliance

!
!##################################################################
!
subroutine tree_resistance(art_resistance,ven_resistance)
!*Descripton:* This subroutine calculates the approximate
! total resistance of a tree so that the solution can be initialised.
! It underestimates the resistance of venous vessels (converging tree)
! as some are added in parallel instead of in a series
    use indices
    use arrays,only: dp,num_elems,elem_field,&
                     elem_cnct,umbilical_inlets,&
                     anastomosis_elem,num_arterial_elems
    use diagnostics, only: enter_exit,get_diagnostics_level
    character(len=60) :: sub_name
!local variables
    real(dp), intent(out) :: art_resistance
    real(dp), intent(out) :: ven_resistance
    real(dp) :: invres,elem_res(num_elems),cap_res
    integer :: num2,ne,ne2,num_connected_elems,inlet_counter,&
               daughter_counter,nu,nc
    integer :: diagnostics_level
    character(len=60) :: vessel_type,rheol_type

    sub_name = 'tree_resistance'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    vessel_type = 'rigid'
    rheol_type = 'constant_visc'


    do nu = 1, num_units
        ne = units(nu)
        nc = elem_cnct(1,1,ne) !capillary unit is downstream of a terminal unit
        call  capillary_resistance(nc,vessel_type, rheol_type,4000.0_dp,3000.0_dp,cap_res,.False.)
    end do

    elem_res(1:num_elems)=elem_field(ne_resist,1:num_elems)

     do ne = 2*num_arterial_elems,num_arterial_elems+1,-1 !Wont work for non matching tree
        invres=0.0_dp
        num_connected_elems = elem_cnct(-1,0,ne) !upstream elems
        if(num_connected_elems.GT.0)then
             daughter_counter = 0
             do num2=1,num_connected_elems
                ne2=elem_cnct(-1,num2,ne)
                if(elem_field(ne_group,ne).eq.2.0_dp)then
                    invres=invres+1.0_dp/elem_res(ne2) !resistance in parallel, for daughter branches
                    daughter_counter = daughter_counter + 1
                end if
             enddo
             if(daughter_counter.GT.0)then
                elem_res(ne)=elem_res(ne)+1.0_dp/invres !resistance in a series
             endif !daughters
          endif !connected
    end do!nw
    ven_resistance =  elem_res(num_arterial_elems+1)

    do ne=num_arterial_elems,1,-1
       invres=0.0_dp
       !exclude the anastomosis elements if ant exists
       if(elem_field(ne_group,ne).ne.3.)then!(anastomosis_elem.EQ.0).OR.(ne.NE.anastomosis_elem))then
          num_connected_elems = elem_cnct(1,0,ne) !downstream elts
          if(num_connected_elems.GT.0)then
             daughter_counter = 0
             do num2=1,num_connected_elems
                ne2=elem_cnct(1,num2,ne)
                if((anastomosis_elem.EQ.0).OR.(ne2.NE.anastomosis_elem))then
                   invres=invres+1.0_dp/elem_res(ne2) !resistance in parallel, for daughter branches
                   daughter_counter = daughter_counter + 1
                endif
             enddo
             if(daughter_counter.GT.0)then
                elem_res(ne)=elem_res(ne)+1.0_dp/invres !resistance in a series
             endif !daughters
          endif !connected
       else
          print *, "excluding anastomosis in total resistance calculation", ne
       endif !not anastomosis
    enddo

    !calculate total tree resistance by summing resistances at each inlet in parallel
    art_resistance = 0
    do inlet_counter=1,count(umbilical_inlets.NE.0)
       art_resistance = art_resistance + 1.0_dp/elem_res(umbilical_inlets(inlet_counter)) !resistance in parallel
    enddo
    art_resistance = 1.0_dp/art_resistance
    call enter_exit(sub_name,2)
end subroutine tree_resistance

end module fetal