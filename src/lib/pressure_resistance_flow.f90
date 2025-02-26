module pressure_resistance_flow
!*Description:* This module contains tools that are used to solve systems of equations representing steady pressure, resistance and flow problems in any branching geometry.
! 
! Descriptions for subroutines that are not included in the subroutine:
!                                                                        
!*calc_sparse_1d_tree:* sets up the system of equations to solve for pressure and flow. 
! It populates the SparseCol,SparseRow, SparceVal and RHS vectors.
!                                                                                 
!*boundary_conditions:* Defines boundary conditions for prq problems
!
  use solve, only: BICGSTAB_LinSolv,pmgmres_ilu_cr
  implicit none
  !Module parameters

  !Module types

  !Module depvar

  !Interfaces
  private
  public evaluate_prq, calculate_stats,tree_resistance,calculate_resistance,capillary_resistance
contains
!###################################################################################
!
subroutine evaluate_prq(mesh_type,bc_type,rheology_type,vessel_type,inlet_flow,inlet_pressure,outlet_pressure)
!*Description:* Solves for pressure and flow in a rigid or compliant tree structure
! Model Types:                                                                     
! mesh_type: can be 'simple_tree' or 'full_plus_tube'. Simple_tree is the input arterial tree 
! without any special features at the terminal level
! 'full_plus_tube' creates a matching venous mesh and has arteries and veins connected by
! capillary units (capillaries are just tubes represented by an element)
!                                                                                  
! boundary condition type (bc_type): pressure or flow
    use indices
    use arrays,only: dp,num_elems,num_nodes,elem_field,elem_nodes,elem_cnct,node_xyz,num_inlets,umbilical_inlets,&
      num_arterial_elems,num_units,units,capillary_model_type
    use diagnostics, only: enter_exit,get_diagnostics_level
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_PRQ" :: EVALUATE_PRQ
  
    character(len=60), intent(in) :: mesh_type,bc_type,rheology_type,vessel_type
    real(dp), intent(in) :: inlet_flow,inlet_pressure,outlet_pressure
    
    !local variables
    integer :: mesh_dof,depvar_types
    integer, allocatable :: mesh_from_depvar(:,:,:)
    integer, allocatable :: depvar_at_node(:,:,:)
    integer, allocatable :: depvar_at_elem(:,:,:)
    integer, dimension(0:2,2) :: depvar_totals
    integer, allocatable :: SparseCol(:)
    integer, allocatable :: SparseRow(:)
    real(dp), allocatable :: SparseVal(:)
    real(dp), allocatable :: RHS(:)
    integer, allocatable :: update_resistance_entries(:)
    integer :: num_vars,NonZeros,MatrixSize
    integer :: AllocateStatus
    integer :: diagnostics_level

    real(dp), allocatable :: prq_solution(:),prq_solution_old(:),solver_solution(:)
    real(dp) :: viscosity,density,inletbc,outletbc,gamma,total_resistance
    logical, allocatable :: FIX(:)
    logical :: ADD=.FALSE.
    character(len=60) :: sub_name
    integer :: no,depvar,nz,ne,SOLVER_FLAG,ne0,ne1,nj,i,iteration_counter,inlet_counter
    real(dp) :: elastance,ERR
    logical :: converged,iterative
    integer :: update_flow_nzz_row
    integer :: np1,np2,nu,nc
    real(dp) :: P1,P2,cap_res
    integer :: inlet_elem
    real(dp) :: inlet_flow_curr,outlet_flow
    character(len=60) :: vessel_type2,rheology_type2

    sub_name = 'evaluate_prq'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

iterative = .false.

if(diagnostics_level.GT.1)then
	print *, "mesh_type=",mesh_type
	print *, "bc_type=",bc_type
	print *, "rheology_type=",rheology_type
	print *, "vessel_type=",vessel_type
endif

if((mesh_type.NE.'full_plus_tube').AND.(mesh_type.NE.'simple_tree'))then
	print *,"unsupported mesh_type",mesh_type
	call exit(1)
endif

if (vessel_type.eq.'rigid') then !Basic model no elasticity
elseif (vessel_type.eq.'elastic') then !Model for vessel elasticity
  elastance=5.0e5_dp!(Pa)
  iterative = .True.
else
   print *,"unsupported vessel_type",vessel_type
   call exit(1)
endif

if (rheology_type.eq.'constant_visc') then !Basic model, no rheology dependent viscosity
elseif (rheology_type.eq.'pries_network') then
   iterative = .True.
elseif(rheology_type.eq.'pries_vessel') then
   iterative = .True.
else
   print *,"unsupported rheology_type",rheology_type
   call exit(1)
endif
elem_field(ne_viscfact,:) = 1.0_dp !initialise viscosity factor

!bc_type: boundary conditions
    !pressure (at inlet and outlets)
    !flow (flow at inlet pressure at outlet).

if((bc_type.NE.'pressure').AND.(bc_type.NE.'flow'))then
	print *,"unsupported bc_type",bc_type
	call exit(1)
elseif((bc_type.EQ.'flow').AND.(inlet_flow.EQ.0))then
	print *, "please set inlet flow"
	call exit(1)
endif



if(bc_type.eq.'pressure')then
    if(inlet_pressure.EQ.0)then
      inletbc=6650.0_dp! Pa (50mmHg) default inlet pressure for human umbilical artery
                       !1 mmHg = 133.322 pascals (Pa)
    else
      inletbc = inlet_pressure
    endif   
elseif(bc_type.eq.'flow')then
    inletbc=inlet_flow
endif

if(outlet_pressure.EQ.0)then
    outletbc=2660.0_dp! Pa (20mmHg) default outlet pressure for human umbilical vein
else
    outletbc = outlet_pressure
endif

if(diagnostics_level.GT.1)then
	print *, "inletbc=",inletbc
	print *, "outletbc",outletbc
endif

if(capillary_model_type.gt.1)then
    iterative = .True.
end if
!!---------PHYSICAL PARAMETERS-----------
density=0.10500e-02_dp !kg/cm3 !density:fluid density
viscosity=0.33600e-02_dp !Pa.s !viscosity: fluid viscosity

!! Allocate memory to depvar arrays
    mesh_dof=num_elems+num_nodes
    depvar_types=2 !pressure/flow
    allocate (mesh_from_depvar(0:2,mesh_dof,0:2), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for mesh_from_depvar array ***"
    allocate (depvar_at_elem(0:2,depvar_types,num_elems), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for depvar_at_elem array ***"
    allocate (depvar_at_node(num_nodes,0:2,depvar_types), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for depvar_at_node array ***"
    allocate (prq_solution(mesh_dof), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for prq_solution array ***"
    prq_solution=0.0_dp !initialise
    if(iterative)then
      allocate (prq_solution_old(mesh_dof), STAT = AllocateStatus)
      if (AllocateStatus /= 0) STOP "*** Not enough memory for prq_solution_old array ***"
    endif
    allocate (update_resistance_entries(num_elems), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for update_resistance array ***"
    update_resistance_entries = 0
    allocate (FIX(mesh_dof), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for FIX array ***"

!! Setting up mappings between nodes, elements and solution depvar
    call calc_depvar_maps(mesh_from_depvar,depvar_at_elem,&
                depvar_totals,depvar_at_node,mesh_dof,num_vars)

!! Define boundary conditions
    !first call to define inlet boundary conditions
    call boundary_conditions(ADD,FIX,bc_type,density,inletbc,outletbc,&
      depvar_at_node,depvar_at_elem,prq_solution,mesh_dof,mesh_type)

    !second call to define pressure bcs at terminal branches
     ADD=.TRUE.
     call boundary_conditions(ADD,FIX,bc_type,density,inletbc,outletbc,&
            depvar_at_node,depvar_at_elem,prq_solution,mesh_dof,mesh_type)

!! Calculate resistance of each element
    call calculate_resistance(viscosity,mesh_type)
    !Initialise capillary resistance (with dummy pressures)
    if(capillary_model_type.gt.1)then
      do nu = 1,num_units
         ne = units(nu)
         nc = elem_cnct(1,1,ne) !capillary unit is downstream of a terminal unit
         call capillary_resistance(nc,vessel_type,rheology_type,4000.0_dp,3000.0_dp,cap_res,.False.)
         elem_field(ne_resist,nc) = cap_res
      enddo
    endif


!! Calculate sparsity structure for solution matrices
    !Determine size of and allocate solution vectors/matrices
    call calc_sparse_size(mesh_dof,FIX,depvar_at_elem,MatrixSize,NonZeros,bc_type)
 	allocate (SparseCol(NonZeros), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for SparseCol array ***"
	allocate (SparseVal(NonZeros), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for SparseVal array ***"
    allocate (SparseRow(MatrixSize+1), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for SparseRow array ***"
    allocate (RHS(MatrixSize), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for RHS array ***"
    allocate (solver_solution(MatrixSize), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    !calculate the sparsity structure
	call calc_sparse_1dtree(bc_type,FIX,mesh_dof,depvar_at_elem, &
        depvar_at_node,NonZeros,MatrixSize,SparseCol,SparseRow,SparseVal,RHS, &
        prq_solution,update_resistance_entries,update_flow_nzz_row)

 !!! Initialise solution vector based on bcs and rigid vessel resistance
   call tree_resistance(total_resistance)
   if(diagnostics_level.GE.2)then
     print *, "total_resistance=",total_resistance
   endif
   call initialise_solution(inletbc,outletbc,(inletbc-outletbc)/total_resistance, &
              mesh_dof,prq_solution,depvar_at_node,depvar_at_elem,FIX)

    !move initialisation to solver solution (skipping BCs).
   no=0
   do depvar=1,mesh_dof !loop over mesh dofs
      if(.NOT.FIX(depvar))then
         no=no+1
         solver_solution(no)=prq_solution(depvar)
      endif
   enddo !mesh_dof
   if((vessel_type.eq."rigid").and.(rheology_type.eq."constant_visc").and.(capillary_model_type.le.1))then
     !! ----CALL SOLVER----
     call pmgmres_ilu_cr(MatrixSize, NonZeros, SparseRow, SparseCol, SparseVal, solver_solution, &
             RHS,  500, 500,1.d-5,1.d-4, SOLVER_FLAG)
     if(SOLVER_FLAG == 0)then
         print *, 'Warning: pmgmres has reached max iterations. Solution may not be valid if this warning persists'
     elseif(SOLVER_FLAG ==2)then
         print *, 'ERROR: pmgmres has failed to converge'
     endif
     !!--TRANSFER SOLVER SOLUTIONS TO FULL SOLUTIONS
     no=0
     do depvar=1,mesh_dof
       if(.NOT.FIX(depvar)) THEN
           no=no+1
		   prq_solution(depvar)=solver_solution(no) !pressure & flow solutions
       endif
     enddo
   else
     iteration_counter = 0
     converged = .False.
     !need an interative loop as problem is non-linear
     do while(.NOT.CONVERGED)
     if (iteration_counter.lt.15)then
       iteration_counter = iteration_counter + 1
       !! ----CALL SOLVER----
       call pmgmres_ilu_cr(MatrixSize, NonZeros, SparseRow, SparseCol, SparseVal, solver_solution, &
              RHS, 500, 500,1.d-5,1.d-4,SOLVER_FLAG)
       if(SOLVER_FLAG == 0)then
         print *, 'Warning: pmgmres has reached max iterations. Solution may not be valid if this warning persists'
       elseif(SOLVER_FLAG ==2)then
          print *, 'ERROR: pmgmres has failed to converge'
       endif
       !!--CALCULATE_ERROR
       ERR=0.0_dp
       no=0
       do depvar=1,mesh_dof
         if(.NOT.FIX(depvar)) THEN
             no=no+1
             prq_solution_old(depvar)=prq_solution(depvar) !temp storage of previous solution
             prq_solution(depvar)=solver_solution(no) !new pressure & flow solutions
             if(DABS(prq_solution(depvar)).GT.0.0_dp)THEN
                ERR=ERR+(prq_solution(depvar)-prq_solution_old(depvar))**2.0_dp/prq_solution(depvar)**2.0_dp
             endif
         endif
       enddo !no2
       ERR=ERR/mesh_dof!sum of error divided by no of unknown depvar

       if(ERR.LE.1.0e-6_dp.AND.(iteration_counter.gt.1))then
         converged=.TRUE.
         print *,"Convergence achieved after",iteration_counter,"iterations",ERR
       else !if error not converged
         print *,"Not converged, error =",ERR
         !!-- UPDATE RESISTANCE OF EACH VESSEL IN THE TREE
         if ((rheology_type.eq.'pries_network').or.(rheology_type.eq.'pries_vessel')) then
            call calculate_viscosity_factor(rheology_type)
         endif
         if(vessel_type.eq.'elastic')then
             !Update vessel radii based on predicted pressures and then update resistance through tree
             call calculate_radius(depvar_at_node,prq_solution,mesh_dof,elastance)
             call calculate_resistance(viscosity,mesh_type)
         endif
       endif !ERR not converged
       if(converged)then
        open(10, file='micro_flow_results.out', status='replace')
       endif
       do nu = 1,num_units
         ne = units(nu)
         nc = elem_cnct(1,1,ne) !capillary unit is downstream of a terminal unit
         np1 = elem_nodes(1,nc)
         np2 = elem_nodes(2,nc)
         P1=prq_solution(depvar_at_node(np1,0,1)) !pressure at start node of capillary element
         P2=prq_solution(depvar_at_node(np2,0,1))!pressure at end node of capillary element
         if(converged)then
           call capillary_resistance(nc,vessel_type,rheology_type,P1,P2,cap_res,.True.)
         else
           call capillary_resistance(nc,vessel_type,rheology_type,P1,P2,cap_res,.False.)
         endif
         elem_field(ne_resist,nc) = cap_res
      enddo
      if(converged)then
        close(10)
      endif

       !!--UPDATE RESISTANCE IN SOLUTION MATRICES
       do ne=1,num_elems !update for all ne
        if(update_resistance_entries(ne).gt.0)then
            nz=update_resistance_entries(ne)
            SparseVal(nz)=-elem_field(ne_resist,ne) !Just updating resistance
        endif
      enddo
      if(bc_type.eq.'flow')then !update RHS to account for element resistance
        do ne=1,num_elems
           depvar = depvar_at_elem(1,1,ne)
           if(FIX(depvar))then
             RHS(update_flow_nzz_row) = prq_solution(depvar)*elem_field(ne_resist,ne)
           endif
        enddo
      endif

     !!--TRANSFER SOLVER SOLUTIONS TO FULL SOLUTIONS
     no=0
     do depvar=1,mesh_dof
        if(.NOT.FIX(depvar)) THEN
            no=no+1
            prq_solution(depvar)=solver_solution(no) !pressure & flow solutions
        endif
     enddo
   endif
  enddo !iterations
 endif




!need to write solution to element/nodal fields for export
    call map_solution_to_mesh(prq_solution,depvar_at_elem,depvar_at_node,mesh_dof)
    !NEED TO UPDATE TERMINAL SOLUTION HERE. LOOP THO' UNITS AND TAKE FLOW AND PRESSURE AT TERMINALS
    call map_flow_to_terminals


     do inlet_counter=1,num_inlets
         inlet_elem = umbilical_inlets(inlet_counter)
         if(inlet_elem.GT.0)then
            inlet_flow_curr = elem_field(ne_Qdot,inlet_elem)
            outlet_flow = elem_field(ne_Qdot,inlet_elem + num_arterial_elems)
            print *, "Flow at inlet element number ",inlet_elem, "(mm3/s) =",inlet_flow_curr
	        print *, "Flow at inlet element number ",inlet_elem, "(ml/min) =",inlet_flow_curr * 0.06_dp
	        print *, "Flow at outlet element number ",inlet_elem + num_arterial_elems, "(mm3/s) =",outlet_flow
	        print *, "Flow at outlet element number ",inlet_elem + num_arterial_elems, "(ml/min) =",outlet_flow * 0.06_dp
	      endif
	 enddo
    deallocate (mesh_from_depvar, STAT = AllocateStatus)
    deallocate (depvar_at_elem, STAT = AllocateStatus)
    deallocate (depvar_at_node, STAT = AllocateStatus)
    deallocate (prq_solution, STAT = AllocateStatus)
    if(iterative)then
      deallocate (prq_solution_old, STAT = AllocateStatus)
    endif
    deallocate (update_resistance_entries,STAT = AllocateStatus)
    deallocate (FIX, STAT = AllocateStatus)
    deallocate (solver_solution, STAT = AllocateStatus)
    deallocate (SparseCol, STAT = AllocateStatus)
    deallocate (SparseVal, STAT = AllocateStatus)
    deallocate (SparseRow, STAT = AllocateStatus)
    deallocate (RHS, STAT = AllocateStatus)
    call enter_exit(sub_name,2)
end subroutine evaluate_prq
!
!###################################################################################
!
subroutine boundary_conditions(ADD,FIX,bc_type,density,inletbc,outletbc,depvar_at_node, &
depvar_at_elem,prq_solution,mesh_dof,mesh_type)

!*Description:* Defines boundary conditions for prq problems

 use arrays,only: dp,num_elems,num_nodes,elem_nodes,elem_cnct,node_xyz,units,&
        num_units
    use diagnostics, only: enter_exit

    integer :: mesh_dof
    integer :: depvar_at_elem(0:2,2,num_elems)
    integer :: depvar_at_node(num_nodes,0:2,2)
    real(dp) :: prq_solution(mesh_dof),inletbc,outletbc,density
    logical:: ADD
    logical :: FIX(mesh_dof)
    character(len=60) ::bc_type,mesh_type

  ! local variables
    integer :: nonode,np,ne,ny1,nj,count_inlet
    character(len=60) :: sub_name

  sub_name = 'boundary_conditions'
  call enter_exit(sub_name,1)

  if(.NOT.ADD)THEN
     ! Initial values
     FIX(1:mesh_dof)=.FALSE.
     prq_solution = 0
     ! Fixed boundary conditions
     ! These are inlet BCs, apply to all inlet BCs (there should only be one)
     count_inlet = 0
     do ne=1,num_elems
        !ne=elems(noelem)
        if (elem_cnct(-1,0,ne) == 0) THEN !Entry element
           count_inlet = count_inlet+1
           if(BC_TYPE == 'pressure')THEN
              np=elem_nodes(1,ne)
              ny1=depvar_at_node(np,1,1) !for fixed pressure BC
              FIX(ny1)=.TRUE. !set fixed
              prq_solution(ny1)=inletbc !Putting BC value into solution array
           else if(BC_TYPE == 'flow')THEN
              ny1=depvar_at_elem(0,1,ne) !fixed
              FIX(ny1)=.TRUE. !set fixed
              prq_solution(ny1)=inletbc !Putting BC value into solution array
           endif
        endif
     enddo
  else !Add terminal pressure BC for all terminal branches
    if(mesh_type.eq.'simple_tree')then
      do nonode=1,num_units
         np=elem_nodes(2,units(nonode)) !Second node in element = terminal node
         ny1=depvar_at_node(np,1,1) !for fixed pressure BC
         FIX(ny1)=.TRUE. !set fixed
         prq_solution(ny1)=outletbc !Putting BC value into solution array
         enddo
     else !if mesh_type is not 'simple_tree'
       do ne=1,num_elems
        !ne=elems(noelem)
        if (elem_cnct(1,0,ne) == 0) THEN !EXIT ELEMENT, no downstream elements exist
              np=elem_nodes(2,ne)
              ny1=depvar_at_node(np,1,1) !for fixed pressure BC
              FIX(ny1)=.TRUE. !set fixed
              prq_solution(ny1)=outletbc !Putting BC value into solution array
          endif
        enddo
     endif
  endif
    call enter_exit(sub_name,2)
  end subroutine boundary_conditions

!
!#################################################################f_###################
!
subroutine calculate_radius(depvar_at_node,prq_solution,mesh_dof,elastance)
    use arrays,only: dp,num_nodes,num_elems,elem_nodes,elem_field
    use indices
    use diagnostics, only: enter_exit,get_diagnostics_level

    integer, intent(in) :: mesh_dof
    integer,intent(in) :: depvar_at_node(num_nodes,0:2,2)
    real(dp),intent(in) ::  prq_solution(mesh_dof)
    real(dp), intent(in) :: elastance
    character(len=60) :: sub_name
    integer :: diagnostics_level

    !local variables
    integer :: ne,nn,np,ny
    real(dp) :: Ptm,h,R0

  sub_name = 'calculate_radius'
  call enter_exit(sub_name,1)
  call get_diagnostics_level(diagnostics_level)
  do ne=1,num_elems
      do nn=1,2
        if(nn.eq.1) np=elem_nodes(1,ne)
        if(nn.eq.2) np=elem_nodes(2,ne)
        ny=depvar_at_node(np,0,1)
        R0 = elem_field(ne_radius,ne)
        if(R0.gt.0.125_dp)then
            h=0.2_dp*R0
        else
            h=0.8_dp*R0
        endif
        Ptm=prq_solution(ny)-11.0_dp*133.0_dp
        if(Ptm.lt.0.0_dp) Ptm=0.0_dp
        if(nn.eq.1) elem_field(ne_radius_in,ne)=R0+3.0_dp*R0**2*Ptm/(4.0_dp*elastance*h)
        if(nn.eq.2) elem_field(ne_radius_out,ne)=R0+3.0_dp*R0**2*Ptm/(4.0_dp*elastance*h)
      enddo

  enddo



  call enter_exit(sub_name,2)
end subroutine calculate_radius

!
!###################################################################################
!
subroutine calculate_resistance(viscosity,mesh_type)
    use arrays,only: dp,num_elems,elem_nodes,node_xyz,&
        elem_field
    use other_consts
    use indices
    use diagnostics, only: enter_exit,get_diagnostics_level
    real(dp):: viscosity
    character(len=60) :: mesh_type
!local variables
    integer :: ne,np1,np2
    real(dp) :: radius,resistance,zeta
    character(len=60) :: sub_name
    integer :: diagnostics_level

  sub_name = 'calculate_resistance'
  call enter_exit(sub_name,1)
  call get_diagnostics_level(diagnostics_level)

!Loop over all elements in model and define resistance for that branch.
    do ne=1,num_elems
       !ne=elems(noelem)
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       radius = (elem_field(ne_radius_in,ne)+elem_field(ne_radius_out,ne))/2.0_dp
       ! element Poiseuille (laminar) resistance in units of Pa.s.mm-3   '
       resistance = 8.0_dp*viscosity*elem_field(ne_viscfact,ne)*elem_field(ne_length,ne)/ &
            (PI*radius**4.0_dp) !laminar resistance

       elem_field(ne_resist,ne) = resistance
       if(diagnostics_level.GT.1)then
       		print *,"TESTING RESISTANCE: element",ne,"resistance",elem_field(ne_resist,ne),"radius",elem_field(ne_radius,ne),&
                    "length",elem_field(ne_length,ne)
       endif

    enddo

    call enter_exit(sub_name,2)
  end subroutine calculate_resistance
!
!##################################################################
!
subroutine calculate_stats(FLOW_GEN_FILE,image_voxel_size,output_level)
!*Descripton:* This subroutine prints the following statistics for the
! feto-placental circulation model:
! arterial, venous, capillary and total vascular volume
! capillary unit surface area
! total resistance of vasculature
! mean, min, max and standard deviation of branch diameter by Strahler order
! mean, min, max and standard deviation of arterial branch diameter by Strahler order
! mean, min, max and standard deviation of venous branch diameter by Strahler order
! coefficient of variation for terminal flow
! terminal flow by generation
!
! Input parameters: FLOW_GEN_FILE - filename to export terminal blood flow statistics
!                                   by generation
!                   image_voxel_size in mm - this is used to calculate the volume
!                                   of vessels in the reconstructed tree that can't be
!                                   resolved in the images

    use indices
    use arrays,only: num_arterial_elems,dp,num_elems, &
                  elem_field,num_units,units,elem_cnct,elem_ordrs, &
                  num_conv,num_conv_gen,cap_resistance,terminal_resistance, &
                  terminal_length, is_capillary_unit,elem_cnct_no_anast, &
                  total_cap_volume,total_cap_surface_area, umbilical_inlets, &
                  umbilical_outlets,elem_nodes,node_field,elem_units_below,num_inlets, &
                  capillary_model_type,num_outlets
    use other_consts, only: PI,MAX_FILENAME_LEN
    use diagnostics, only: enter_exit,get_diagnostics_level
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CALCULATE_STATS" :: CALCULATE_STATS

    character(len=MAX_FILENAME_LEN), intent(in) :: FLOW_GEN_FILE
    real(dp), intent(in) :: image_voxel_size
    integer, intent(in) :: output_level
  !local variables
    integer :: ne,nu,nc,order,max_strahler,no_branches, &
               ne_order,branch,ven_elems,max_gen,np1,np2, &
               inlet_counter,inlet_elem,outlet_elem,n,num_units1,&
               outlet_counter
    real(dp) :: total_vasc_volume, &
                mean_diameter,std_diameter,total_resistance,std_terminal_flow, &
                cof_var_terminal_flow,mean_terminal_flow,small_vessel_volume, diameter, &
                cap_length,cof_var_cap_flow,inlet_flow,inlet_pressure,outlet_pressure, &
                resistance,volume_fed_by_inlet,mean_terminal_flow1,mean_terminal_flow2,&
                std_terminal_flow1,std_terminal_flow2,total_surface_area
    real(dp) :: arterial_vasc_volume,arterial_surface_area,arteriole_vasc_volume, &
         arteriole_surface_area,cap_vasc_volume,cap_surface_area,venule_vasc_volume,&
         venule_surface_area,vein_vasc_volume,vein_surface_area,radius,sum_cap_resistance,&
         sum_terminal_lengths
    real(dp) :: outlet_flow
    integer :: strahler_orders(num_elems)
    integer :: generations(num_elems)
    real(dp),allocatable :: diameter_by_strahler(:,:)
    real(dp),allocatable :: art_diameter_by_strahler(:,:)
    real(dp),allocatable :: ven_diameter_by_strahler(:,:)
    real(dp),allocatable :: terminal_flow_by_gen(:,:)
    integer,allocatable :: branch_count(:)
    integer,allocatable :: gen_branch_count(:)
    logical :: elems_under_inlet1(num_arterial_elems)
    logical :: elems_under_inlet2(num_arterial_elems)
    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'calculate_stats'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

  !Output level == 2 then we have no capillaries
  !calculate arterial vascular volume
   arterial_vasc_volume = 0.0_dp
   arterial_surface_area = 0.0_dp
   if (output_level.eq.1)then
       arteriole_vasc_volume = 0.0_dp
       arteriole_surface_area = 0.0_dp
       cap_vasc_volume = 0.0_dp
       cap_surface_area = 0.0_dp
       venule_vasc_volume = 0.0_dp
       venule_surface_area = 0.0_dp
       vein_vasc_volume = 0.0_dp
       vein_surface_area = 0.0_dp
       sum_cap_resistance = 0.0_dp
       sum_terminal_lengths = 0.0_dp
       mean_terminal_flow = 0.0_dp
   endif
   do ne=1,num_elems
       !RADIUS IS AVERAGE OF THE TWO RADII
       radius = (elem_field(ne_radius_in,ne)+elem_field(ne_radius_out,ne))/2.0_dp
      if(elem_field(ne_group,ne).eq.0)then
        arterial_vasc_volume = arterial_vasc_volume + PI*radius**2.0_dp*elem_field(ne_length, ne)
        arterial_surface_area = arterial_surface_area  + PI* 2.0_dp*radius*elem_field(ne_length, ne)
      else if((elem_field(ne_group,ne).eq.1).and.(output_level.eq.1))then
         arteriole_vasc_volume = arteriole_vasc_volume + elem_field(ne_artvol,ne)
      	 arteriole_surface_area = arteriole_surface_area  + elem_field(ne_artsa,ne)
         cap_vasc_volume = cap_vasc_volume + elem_field(ne_vol,ne)
      	 cap_surface_area = cap_surface_area  + elem_field(ne_sa,ne)
         venule_vasc_volume = venule_vasc_volume + elem_field(ne_veinvol,ne)
      	 venule_surface_area = venule_surface_area  + elem_field(ne_veinsa,ne)
         if(capillary_model_type.gt.1)then
            sum_cap_resistance = sum_cap_resistance + elem_field(ne_resist,ne)
            sum_terminal_lengths = sum_terminal_lengths + elem_field(ne_length,ne)
         end if
         mean_terminal_flow = mean_terminal_flow + elem_field(ne_Qdot,ne)

      else if((elem_field(ne_group,ne).eq.2).and.(output_level.eq.1))then
        vein_vasc_volume = vein_vasc_volume + PI*radius**2.0_dp*elem_field(ne_length, ne)
        vein_surface_area = vein_surface_area  + PI* 2.0_dp*radius*elem_field(ne_length, ne)
      else
        !Anastomosis
      endif
   enddo
   mean_terminal_flow = mean_terminal_flow/dble(num_units)
   print *, "##########################    ARTERIES   #########################"
   print *, "Arterial elems count",num_arterial_elems
   print *, "Arterial vascular volume (cm**3/ml) = ",arterial_vasc_volume/1000.0_dp !mm3 to cm3
   print*, "Arterial surface area (cm**2) = ", arterial_surface_area/100.0_dp
   print*, "Arterial surface area (m**2) = ", arterial_surface_area/(100.0_dp*10000.0_dp)
   print *, "##########################    ARTERIES   #########################"
   if(output_level.eq.1) then
       print *, "###########################  ARTERIOLE  ###########################"
       print *, "Arteriole vascular volume (cm**3/ml) = ",arteriole_vasc_volume/1000.0_dp !mm3 to cm3
       print*, "Arteriole surface area (cm**3) = ", arteriole_surface_area/100.0_dp
       print*, "Arteriole surface area (m**3) = ", arteriole_surface_area/(100.0_dp*10000.0_dp)
       print *, "###########################  ARTERIOLE  ###########################"
       print *, "##########################   CAPILLARY   ##########################"
       print *, "Capillary elems count",num_units
       print *, "Capillary vascular volume (cm**3/ml) = ",cap_vasc_volume/1000.0_dp !mm3 to cm3
       print*, "Capillary surface area (cm**3) = ", cap_surface_area/100.0_dp
       print*, "Capillary surface area (m**3) = ", cap_surface_area/(100.0_dp*10000.0_dp)
       print *, "##########################   CAPILLARY   ##########################"
       print *, "############################  VENULE  #############################"
       print *, "Venule vascular volume (cm**3/ml) = ",venule_vasc_volume/1000.0_dp !mm3 to cm3
       print*, "Venule surface area (cm**3) = ", venule_surface_area/100.0_dp
       print*, "Venule surface area (m**3) = ", venule_surface_area/(100.0_dp*10000.0_dp)
       print *, "############################  VENULE  #############################"
       print *, "###########################   VEINS  #############################"
       print *, "Vein vascular volume (cm**3/ml) = ",vein_vasc_volume/1000.0_dp !mm3 to cm3
       print*, "Vein surface area (cm**3) = ", vein_surface_area/100.0_dp
       print*, "Vein surface area (m**3) = ", vein_surface_area/(100.0_dp*10000.0_dp)
       print *, "############################   VEINS  #############################"
   endif

   !Print capillary convolute number, generations and total length of capillaries
   if(capillary_model_type.eq.1)then
      print *, "Resistance of individual capillary conduits (Pa.s/mm**3) =",cap_resistance
      print *, "Resistance of terminal unit (Pa.s/mm**3) =",terminal_resistance
   else
      print *, "Baseline resistance of individual capillary conduits (Pa.s/mm**3) =",cap_resistance
      print *, "Average (mean) resistance of terminal unit (Pa.s/mm**3) =", sum_cap_resistance/dble(num_units)

   end if


   !calculate total vascular resistance (Pressure in - Pressure out)/Blood Flow in
   inlet_flow = 0.0_dp
   inlet_pressure = 0.0_dp

   do inlet_counter=1,num_inlets
      inlet_elem = umbilical_inlets(inlet_counter)
      if(inlet_elem.GT.0)then
         inlet_flow = inlet_flow + elem_field(ne_Qdot,inlet_elem)
         !get the first node
         np1 = elem_nodes(1,inlet_elem)
         inlet_pressure = node_field(nj_bv_press,np1)
      endif
   enddo
   outlet_pressure = 0.0_dp
   if(allocated(umbilical_outlets))then
      outlet_elem = umbilical_outlets(1)
   else !get the first unit
      outlet_elem = units(1)
   endif
   !get the second node
   np2 = elem_nodes(2,outlet_elem)
   outlet_pressure = node_field(nj_bv_press,np2)
   resistance = (inlet_pressure - outlet_pressure)/inlet_flow
   print *, "######################## OVERALL PRESSURE/FLOW #####################"
   print *, "Inlet pressure (Pa) =", inlet_pressure
   print *, "Inlet pressure (mmHg) = ", inlet_pressure/133.322_dp
   print *, "Outlet pressure (Pa) =", outlet_pressure
   print *, "Outlet pressure (mmHg) =", outlet_pressure/133.322_dp
   print *, "Flow (sum of all inlet flows) (mm3/s) =", inlet_flow
   print *, "Flow (sum of all inlet flows) (ml/min) =", inlet_flow * 0.06_dp
   print *, "Total vascular resistance (Pin-Pout)/Flow (Pa.s/mm**3) = ",resistance
   print *, "######################## OVERALL PRESSURE/FLOW #####################"

   !mean, min, max and std of branch diameter by Strahler order
   if(output_level.gt.0)then
     strahler_orders(:) = elem_ordrs(no_sord, :)
     max_strahler = maxval(strahler_orders)
     allocate(diameter_by_strahler(max_strahler,num_elems))
     allocate(branch_count(max_strahler))
     diameter_by_strahler = 0
     branch_count = 0

     do ne=1,num_elems
        ne_order = strahler_orders(ne)
        if(ne_order.ge.1)then
           no_branches = branch_count(ne_order)
           no_branches = no_branches + 1
           diameter_by_strahler(ne_order,no_branches) = elem_field(ne_radius,ne) * 2
           branch_count(ne_order) = no_branches
       endif
     enddo
     print *, "Vessel diameter (mm) by Strahler order:"
     print *, "Strahler_order,number_of_elements,mean_diameter,min_diameter,max_diameter,std"
     do order=1, max_strahler
        mean_diameter = sum(diameter_by_strahler(order,:))/branch_count(order)
        std_diameter = 0
        do branch=1,branch_count(order)
          std_diameter = std_diameter + (diameter_by_strahler(order,branch) - mean_diameter)**2
        enddo
        std_diameter = std_diameter/branch_count(order)
        std_diameter = SQRT(std_diameter)
        print *, order,",",branch_count(order),",",branch_count(order),",",mean_diameter,",", &
            minval(diameter_by_strahler(order,:),MASK = diameter_by_strahler(order,:) .GT.0),",", &
            maxval(diameter_by_strahler(order,:)),",",std_diameter
     enddo
     if(output_level.eq.1)then
         !Arterial vessel diameter by Strahler order
         allocate(art_diameter_by_strahler(max_strahler,num_arterial_elems))
         art_diameter_by_strahler = 0
         branch_count = 0
         do ne=1,num_arterial_elems
            ne_order = strahler_orders(ne)
            if(ne_order.ge.1)then
             no_branches = branch_count(ne_order)
             no_branches = no_branches + 1
             art_diameter_by_strahler(ne_order,no_branches) = elem_field(ne_radius,ne) * 2
             branch_count(ne_order) = no_branches
            endif
         enddo
         print *, "Arterial vessel diameter (mm) by Strahler order:"
         print *, "Strahler_order,number_of_elements,mean_diameter,min_diameter,max_diameter,std"
         do order=1, max_strahler
            mean_diameter = sum(art_diameter_by_strahler(order,:))/branch_count(order)
            std_diameter = 0
            do branch=1,branch_count(order)
              std_diameter = std_diameter + (art_diameter_by_strahler(order,branch) - mean_diameter)**2
            enddo
            std_diameter = std_diameter/branch_count(order)
            std_diameter = SQRT(std_diameter)
            print *, order,",",branch_count(order),",",mean_diameter,",", &
               minval(art_diameter_by_strahler(order,:),MASK = art_diameter_by_strahler(order,:) .GT.0),",", &
                maxval(art_diameter_by_strahler(order,:)),",",std_diameter
         enddo

         !Venous vessel diameter by Strahler order
         !all elements - arterial elements - number of capillaries (the same as number of terminal units)
         ven_elems = num_elems-num_arterial_elems-num_units
         write(*,*) ven_elems,num_elems,num_arterial_elems,num_units,max_strahler
         if (ven_elems.GT.0)then!checking there are venous elements
          allocate(ven_diameter_by_strahler(max_strahler,ven_elems))
          ven_diameter_by_strahler = 0
          branch_count = 0
          do ne=num_arterial_elems+1,num_elems
             if(elem_field(ne_group,ne).eq.2.0_dp)then !vein
                ne_order = strahler_orders(ne)
	        if(ne_order.ge.1)then
                   no_branches = branch_count(ne_order)
                   no_branches = no_branches + 1
                   ven_diameter_by_strahler(ne_order,no_branches) = elem_field(ne_radius,ne) * 2
                   branch_count(ne_order) = no_branches
	        endif
             endif
          enddo
          print *, "Venous vessel diameter (mm) by Strahler order:"
          print *, "Strahler_order,number_of_elements,mean_diameter,min_diameter,max_diameter,std"
          do order=1, max_strahler
             mean_diameter = sum(ven_diameter_by_strahler(order,:))/branch_count(order)
             std_diameter = 0
             do branch=1,branch_count(order)
                std_diameter = std_diameter + (ven_diameter_by_strahler(order,branch) - mean_diameter)**2
             enddo
             std_diameter = std_diameter/branch_count(order)
             std_diameter = SQRT(std_diameter)
             print *, order,",",branch_count(order),",",mean_diameter,",", &
               minval(ven_diameter_by_strahler(order,:),MASK = ven_diameter_by_strahler(order,:) .GT.0),",", &
                maxval(ven_diameter_by_strahler(order,:)),",",std_diameter
          enddo
         endif !(ven_elems.GT.0)
       endif!outputlevel
     end if

  !coefficient of variation for terminal flow
   if(output_level.eq.2)then
       mean_terminal_flow = 0.0_dp
       do nu=1,num_units
           ne = units(nu)
           mean_terminal_flow = mean_terminal_flow + elem_field(ne_Qdot,ne)
       enddo
       mean_terminal_flow= mean_terminal_flow/num_units

   end if
   std_terminal_flow = 0
   do nu=1,num_units
      ne = units(nu)
      std_terminal_flow = std_terminal_flow + (elem_field(ne_Qdot,ne) - mean_terminal_flow)**2
   enddo
   std_terminal_flow = std_terminal_flow/dble(num_units)
   std_terminal_flow = SQRT(std_terminal_flow)
   cof_var_terminal_flow = std_terminal_flow/mean_terminal_flow
   print *, "Coefficient of variation for terminal flow (%) = ", cof_var_terminal_flow * 100.0_dp
   print *, "Mean terminal flow (mm3/s) = ",mean_terminal_flow
   print *, "Standard deviation of terminal flow (mm3/s) = ",std_terminal_flow

   !terminal flow by generation
   generations(:) = elem_ordrs(no_gen, :)
   max_gen = maxval(generations)
   allocate(terminal_flow_by_gen(max_gen,num_units))
   allocate(gen_branch_count(max_gen))
   terminal_flow_by_gen = 0
   gen_branch_count = 0

   !print all terminal flows and their corresponding generations to a file
   open(10, file=FLOW_GEN_FILE, status="replace")
   write(10,*) 'terminal_blood_flow,generation'
   do nu=1,num_units
      ne = units(nu)
      ne_order = generations(ne)
      no_branches = gen_branch_count(ne_order)
      no_branches = no_branches + 1
      terminal_flow_by_gen(ne_order,no_branches) = elem_field(ne_Qdot,ne)
      gen_branch_count(ne_order) = no_branches
      write(10,*) elem_field(ne_Qdot,ne),',',ne_order
   enddo
   close(10)
   if(output_level.gt.0) then
     print *, "Terminal flow (mm**3/s) by generation:"
     print *, "Generation,number_of_terminal_units,mean_flow,min_flow,max_flow,std"
     do order=1, max_gen
      if(gen_branch_count(order).GT.0)then
         mean_terminal_flow = sum(terminal_flow_by_gen(order,:))/gen_branch_count(order)
         std_terminal_flow = 0
         do branch=1,gen_branch_count(order)
            std_terminal_flow = std_terminal_flow + (terminal_flow_by_gen(order,branch) - mean_terminal_flow)**2
         enddo
         std_terminal_flow = std_terminal_flow/gen_branch_count(order)
         std_terminal_flow = SQRT(std_terminal_flow)
         print *, order,",",gen_branch_count(order),",",mean_terminal_flow,",", &
            minval(terminal_flow_by_gen(order,:),MASK = terminal_flow_by_gen(order,:) .GT.0),",", &
            maxval(terminal_flow_by_gen(order,:)),",",std_terminal_flow
      else
	     print *,order,",0,0,0,0,0"
      endif
     enddo
   endif

   !if more than one inlet
   if(count(umbilical_inlets.NE.0).GT.1)then

      !print pressure and flow at each inlet
      do inlet_counter=1,num_inlets
         inlet_elem = umbilical_inlets(inlet_counter)
         if(inlet_elem.GT.0)then
            inlet_flow = elem_field(ne_Qdot,inlet_elem)
            print *, "Flow at inlet element number ",inlet_elem, "(mm3/s) =",inlet_flow
	        print *, "Flow at inlet element number ",inlet_elem, "(ml/min) =",inlet_flow * 0.06_dp

            !get the first node
            np1 = elem_nodes(1,inlet_elem)
            inlet_pressure = node_field(nj_bv_press,np1)
  	    print *, "Pressure at inlet node",np1,"inlet element",inlet_elem,"(Pa) =", inlet_pressure
	    print *, "Pressure at inlet node",np1,"inlet element",inlet_elem,"(mmHg) =", inlet_pressure/133.322_dp
         endif
      enddo
     do outlet_counter=1,num_outlets
          outlet_elem = umbilical_outlets(outlet_counter)
         if(outlet_elem.GT.0)then
            outlet_flow = elem_field(ne_Qdot,outlet_elem)
            print *, "Flow at outlet element number ",outlet_elem, "(mm3/s) =",outlet_flow
	        print *, "Flow at outlet element number ",outlet_elem, "(ml/min) =",outlet_flow * 0.06_dp

            !get the first node
            np1 = elem_nodes(1,outlet_elem)
            outlet_pressure = node_field(nj_bv_press,np1)
  	        print *, "Pressure at outlet node",np1,"outlet element",outlet_elem,"(Pa) =", outlet_pressure
	        print *, "Pressure at outlet node",np1,"outlet element",outlet_elem,"(mmHg) =", outlet_pressure/133.322_dp
         end if
      end do
      do inlet_counter=1,num_inlets
         inlet_elem = umbilical_inlets(inlet_counter)
	 print *, "Number of terminal units below inlet element",inlet_elem, "=",elem_units_below(inlet_elem)
      enddo

      elems_under_inlet1 = .FALSE.
      elems_under_inlet2 = .FALSE.
      !populate arrays of arterial elements under each inlet
      !elements directly downstream of the inlets
      do inlet_counter=1,num_inlets
         inlet_elem = umbilical_inlets(inlet_counter)
	 do n=1,elem_cnct_no_anast(1,0,inlet_elem) !Assumes two inlets
            if(inlet_counter.EQ.1)then
               elems_under_inlet1(elem_cnct_no_anast(1,n,inlet_elem)) = .TRUE.
            else
               elems_under_inlet2(elem_cnct_no_anast(1,n,inlet_elem)) = .TRUE.
            endif
         enddo
      enddo
      do ne=1,num_arterial_elems
         if(ALL(umbilical_inlets.NE.ne))then
            !check which inlet the upstream elements are under and assign this element to the same inlet
            do n=1,elem_cnct_no_anast(-1,0,ne)
               if(elems_under_inlet1(elem_cnct_no_anast(-1,n,ne)))then
	          elems_under_inlet1(ne) = .TRUE.
               elseif(elems_under_inlet2(elem_cnct_no_anast(-1,n,ne)))then
                  elems_under_inlet2(ne) = .TRUE.
               endif
            enddo
         endif
      enddo
      !print arterial volume fed by each inlet
      !volume under the first inlet
      volume_fed_by_inlet = 0.0_dp
      do ne=1,num_arterial_elems
         if(elems_under_inlet1(ne))then
            volume_fed_by_inlet = volume_fed_by_inlet + elem_field(ne_vol,ne)
         endif
      enddo
      print *,"Arterial volume fed by inlet element",umbilical_inlets(1),"(cm3) =",volume_fed_by_inlet/1000
      volume_fed_by_inlet = 0.0_dp
      do ne=1,num_arterial_elems
         if(elems_under_inlet2(ne))then
            volume_fed_by_inlet = volume_fed_by_inlet + elem_field(ne_vol,ne)
         endif
      enddo
      print *,"Arterial volume fed by inlet element",umbilical_inlets(2),"(cm3) =",volume_fed_by_inlet/1000

      !coefficient of variation of terminal flow under each inlet
      !standard deviation and mean flow under each inlet
      mean_terminal_flow1 = 0
      mean_terminal_flow2 = 0
      num_units1 = 0
      do nu=1,num_units
         ne = units(nu)
         if(elems_under_inlet1(ne))then
            mean_terminal_flow1 = mean_terminal_flow1 + elem_field(ne_Qdot,ne)
            num_units1 = num_units1 + 1
         elseif(elems_under_inlet2(ne))then
            mean_terminal_flow2 = mean_terminal_flow2 + elem_field(ne_Qdot,ne)
         endif
      enddo
      mean_terminal_flow1 = mean_terminal_flow1/dble(num_units1)
      mean_terminal_flow2 = mean_terminal_flow2/(num_units - num_units1)

      std_terminal_flow1 = 0
      std_terminal_flow2 = 0
      do nu=1,num_units
         ne = units(nu)
         if(elems_under_inlet1(ne))then
            std_terminal_flow1 = std_terminal_flow1 + (elem_field(ne_Qdot,ne) - mean_terminal_flow1)**2
         elseif(elems_under_inlet2(ne))then
            std_terminal_flow2 = std_terminal_flow2 + (elem_field(ne_Qdot,ne) - mean_terminal_flow2)**2
         endif
      enddo
      std_terminal_flow1 = std_terminal_flow1/dble(num_units1)
      std_terminal_flow1 = SQRT(std_terminal_flow1)
      cof_var_terminal_flow = std_terminal_flow1/mean_terminal_flow1
      print *, "Coefficient of variation for terminal flow (%) under inlet element", &
                              umbilical_inlets(1),"=",cof_var_terminal_flow * 100
      print *, "Mean terminal flow (mm3/s) under inlet element",umbilical_inlets(1),"=",mean_terminal_flow1
      print *, "Standard deviation of terminal flow (mm3/s) under inlet element",umbilical_inlets(1),"=",&
              std_terminal_flow1

      std_terminal_flow2 = std_terminal_flow2/(num_units - num_units1)
      std_terminal_flow2 = SQRT(std_terminal_flow2)
      cof_var_terminal_flow = std_terminal_flow2/mean_terminal_flow2
      print *, "Coefficient of variation for terminal flow (%) under inlet element", &
                              umbilical_inlets(2),"=",cof_var_terminal_flow * 100
      print *, "Mean terminal flow (mm3/s) under inlet element",umbilical_inlets(2),"=",mean_terminal_flow2
      print *, "Standard deviation of terminal flow (mm3/s) under inlet element",umbilical_inlets(2),"=",&
                std_terminal_flow2

   endif

   call enter_exit(sub_name,2)

end subroutine calculate_stats


!
!###################################################################################
!
subroutine calculate_viscosity_factor(rheology_type)
    use arrays
    use other_consts
    use indices
    use diagnostics, only: enter_exit,get_diagnostics_level

    character(len=60), intent(in) :: rheology_type

    !local variables
    integer :: ne,nc,ne_d1,ne_d22
    real(dp) :: controlhb,beta,radius

    character(len=60) :: sub_name
    integer :: diagnostics_level
    sub_name = 'calculate_viscosity_factor'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)


    elem_field(ne_hb,:) = 0.45_dp !Normal_fetal hematocrit
    !controlhb = 0.45_dp

    do ne = 1,num_elems
      radius = (elem_field(ne_radius_in,ne) + elem_field(ne_radius_out,ne))/2.0_dp*1000.0_dp !mm to um
      call viscosity_from_radius(radius,elem_field(ne_hb,ne),elem_field(ne_viscfact,ne))
    enddo

    call enter_exit(sub_name,2)
  end subroutine calculate_viscosity_factor
!
!##################################################################
!
subroutine calc_depvar_maps(mesh_from_depvar,depvar_at_elem,depvar_totals,depvar_at_node,mesh_dof,num_vars)
!*Description:* This subroutine calculates the mapping between the nodes and elements and
! the problem dependent variables that are needed for matrix setup and solution

    use arrays,only: num_elems,num_nodes,elem_nodes
    use diagnostics, only: enter_exit,get_diagnostics_level
    character(len=60) :: sub_name

     integer :: mesh_dof
     integer :: mesh_from_depvar(0:2,mesh_dof,0:2)
     integer :: depvar_at_elem(0:2,2,num_elems)
     integer :: depvar_at_node(num_nodes,0:2,2)
     integer :: depvar_totals(0:2,2)
     integer :: num_vars
!     local variables
    integer :: ny_start=0
    integer :: nc,ne,nn,np,nrc,ny
    integer :: diagnostics_level

    sub_name = 'calc_depvar_maps'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

     depvar_totals = 0
     mesh_from_depvar = 0
     depvar_at_elem = 0
     depvar_at_node = 0

!nrc = loops from 0,1,2

!  Set up mapping arrays for current region:
!  includes nested loop creating depvar_at_node & depvar_at_elem consecutively depvar_at_node(np,nrc,nc) depvar_at_elem(nrc,nc,ne)

     do nrc=0,2 !row or column no
        ny=0 !depvar tag
        do nc=1,2 !no of dependent depvar types
          ! if(nrc.NE.0) ny=ny_start !--> resets to ny_start only for nrc=1,2. Maybe should also do for nrc=1??
           ny=ny_start
           do ne=1,num_elems
              if((diagnostics_level.GT.1).AND.(nc.NE.2).AND.(nrc.NE.2)) THEN
              	print *,""
              	print *,"element ne",ne
              endif
              !ne=elems(noelem)
              do nn=1,2 !Loop over number of nodes in each element: always=2 in 1D element
                 np=elem_nodes(nn,ne)
                    if(depvar_at_node(np,nrc,nc).EQ.0)THEN
                    !Checks if this node already been done
                    ny=ny+1
                     if(ny.GT.mesh_dof) THEN
                       print *,"1. Need to increase mesh_dof! ny=",ny
                       call exit(1)
                    endif
                    if(nrc.NE.0) THEN
                       if(ny.GT.depvar_totals(nrc,nc)) depvar_totals(nrc,nc)=ny
                    endif
                    depvar_at_node(np,nrc,nc)=ny
                    if(nrc.NE.1.OR.nc.EQ.1) THEN
                       ! don't set mesh_from_depvar for nrc=1(rows), nc<>1(LHS) cases
                       mesh_from_depvar(0,ny,nrc)=1 !mesh dof is nodal based
                       mesh_from_depvar(1,ny,nrc)=np
                       mesh_from_depvar(2,ny,nrc)=nc
                    endif !nrc.NE.1
                 endif !depvar_at_node.EQ.0

                if((diagnostics_level.GT.1).AND.(nc.NE.2).AND.(nrc.NE.2)) THEN !arrays where nc=2 or nrc=2 don't seem to be used anywhere
                		print *,"depvar_at_node(",np,",",nrc,",",nc,")=", depvar_at_node(np,nrc,nc)
                endif
              enddo !nn (np)
              ny=ny+1
              if(ny.GT.mesh_dof) THEN
                 print *,"2. Need to increase mesh_dof! ny=",ny
                 call exit(1)
              endif
              if(nrc.NE.0) THEN
                 if(ny.GT.depvar_totals(nrc,nc)) depvar_totals(nrc,nc)=ny
              endif
              depvar_at_elem(nrc,nc,ne)=ny
              if(nrc.NE.1.OR.nc.EQ.1) THEN
                 !                 don't set mesh_from_depvar for nrc=1(rows), nc<>1(LHS) cases
                 mesh_from_depvar(0,ny,nrc)=2 !mesh dof is element based
                 mesh_from_depvar(1,ny,nrc)=nc
                 mesh_from_depvar(2,ny,nrc)=ne
              endif

              if((diagnostics_level.GT.1).AND.(nc.NE.2).AND.(nrc.NE.2)) THEN !arrays where nc=2 or nrc=2 don't seem to be used anywhere
              	 print *,"depvar_at_elem(",nrc,",",nc,",",ne,")=", depvar_at_elem(nrc,nc,ne)
              endif
           enddo !noelem (ne)

        enddo !nc
     enddo !nrc
     num_vars=ny
     if(diagnostics_level.GT.1)then
     	print *,"Max ny number (number of variables)=",ny
	 endif
    call enter_exit(sub_name,2)
end subroutine calc_depvar_maps
!
!##################################################################
!
subroutine tree_resistance(resistance)
!*Descripton:* This subroutine calculates the approximate
! total resistance of a tree so that the solution can be initialised.
! It underestimates the resistance of venous vessels (converging tree)
! as some are added in parallel instead of in a series
    use indices
    use arrays,only: dp,num_elems,elem_field,&
                     elem_cnct,umbilical_inlets,&
                     anastomosis_elem
    use diagnostics, only: enter_exit,get_diagnostics_level
    character(len=60) :: sub_name
!local variables
    real(dp), intent(out) :: resistance
    real(dp) :: invres,elem_res(num_elems)
    integer :: num2,ne,ne2,num_connected_elems,inlet_counter,&
               daughter_counter
    integer :: diagnostics_level

    sub_name = 'tree_resistance'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    elem_res(1:num_elems)=elem_field(ne_resist,1:num_elems)

    do ne=num_elems,1,-1
       invres=0.0_dp
       !exclude the anastomosis elements if ant exists
       if(elem_field(ne_group,ne).ne.3.)then!(anastomosis_elem.EQ.0).OR.(ne.NE.anastomosis_elem))then
          num_connected_elems = elem_cnct(1,0,ne)
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
    resistance = 0
    do inlet_counter=1,count(umbilical_inlets.NE.0)
       resistance = resistance + 1.0_dp/elem_res(umbilical_inlets(inlet_counter)) !resistance in parallel
    enddo
    resistance = 1.0_dp/resistance
    if(diagnostics_level.GT.0)then
       print *,"tree resistance to initialise solution: ",resistance
    endif

    call enter_exit(sub_name,2)
end subroutine tree_resistance
!
!##################################################################
!
subroutine initialise_solution(pressure_in,pressure_out,cardiac_output,mesh_dof,prq_solution,depvar_at_node,depvar_at_elem,FIX)
!*Description:* This subroutine calculates an estimate for initial solution to a prq problem
! based on cardiact output and pressure BCs. Unknown pressure variables are set to the average
! of incoming and outgoing pressures.  Unknown flow variables are set to the cardiac output
! divided by 2 to the power of number of branching generations

    use indices
    use arrays,only: dp,num_elems,elem_ordrs,elem_nodes,num_nodes
    use diagnostics, only: enter_exit
    integer, intent(in) :: mesh_dof
    integer,intent(in) :: depvar_at_elem(0:2,2,num_elems)
    integer,intent(in) :: depvar_at_node(num_nodes,0:2,2)
    real(dp), intent(in) :: pressure_in, pressure_out,cardiac_output
    real(dp) :: prq_solution(mesh_dof)
    logical, intent(in) :: FIX(mesh_dof)
!local variables
    integer :: nn,ne,np,n_depvar
    character(len=60) :: sub_name
   sub_name = 'initialise_solution'
    call enter_exit(sub_name,1)
    do ne=1,num_elems
       !ne=elems(noelem)
         do nn=1,2 !Loop over number of nodes in each element: always=2 in 1D element
         np=elem_nodes(nn,ne) !Node number
         n_depvar=depvar_at_node(np,0,1) !--> This will be a pressure because it's the depvar at a node
         if(.NOT.FIX(n_depvar))then
           prq_solution(n_depvar)=(pressure_in+pressure_out)/2.0
         endif
         n_depvar=depvar_at_elem(0,1,ne) !--> This will be a flow because it's the depvar at the element
         if(.NOT.FIX(n_depvar))then
           prq_solution(n_depvar)=cardiac_output/(2.0_dp**(elem_ordrs(1,ne)-1)) !Here you can use the generation number to split flow
         endif
       enddo
    enddo
    call enter_exit(sub_name,2)
  end subroutine initialise_solution
!
!##################################################################
!
subroutine calc_sparse_size(mesh_dof,FIX,depvar_at_elem,MatrixSize,NonZeros,bc_type)
!*Description:* This subroutine calculates sparsity sizes

    use diagnostics, only: enter_exit,get_diagnostics_level
    use arrays,only: num_elems,elem_nodes,num_nodes,elems_at_node

    integer, intent(in) :: mesh_dof
    logical, intent(in) :: FIX(mesh_dof)
    integer,intent(in) :: depvar_at_elem(0:2,2,num_elems)
    integer,intent(inout) :: MatrixSize
    integer,intent(inout) :: NonZeros
    character(len=60), intent(in) :: bc_type
!local variables
    integer :: i,ne,np,fixed_variables, fixed_flows, fixed_pressures
    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'calc_sparse_size'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

 	fixed_variables = 0
 	!count fixed variables
 	do i=1,mesh_dof
 		if(FIX(i))then
 			fixed_variables = fixed_variables + 1
 		endif
 	enddo
 	MatrixSize = mesh_dof - fixed_variables

 	!get count of fixed flows
 	fixed_flows = 0
 	do ne=1,num_elems
 		if(FIX(depvar_at_elem(1,1,ne)))then
 			fixed_flows = fixed_flows + 1
 		endif
 	enddo

 	fixed_pressures = fixed_variables - fixed_flows

 	!count of pressure equations = (number of elements * 3 variables in each equation) - fixed pressures - fixed flows
 	NonZeros = num_elems*3 - fixed_pressures - fixed_flows
 	!count of conservation of flow equations = sum of elements connected to nodes which have at least 2 connected elements - fixed flows
 	do np=1, num_nodes
 		if(elems_at_node(np,0).GT.1)then
 			NonZeros = NonZeros + elems_at_node(np,0)
 		endif
 	enddo
 	NonZeros = NonZeros - fixed_flows


 	if(diagnostics_level.GT.1)then
 		print *,"MatrixSize",MatrixSize
 		print *,"NonZeros",NonZeros
 	endif

end subroutine calc_sparse_size

!##################################################################
!
subroutine map_solution_to_mesh(prq_solution,depvar_at_elem,depvar_at_node,mesh_dof)
!*Description:* This subroutine maps the solution array to appropriate nodal and element fields

    use indices
    use arrays,only: dp,num_nodes,num_elems,elem_field,node_field
    use diagnostics, only: enter_exit,get_diagnostics_level

    integer, intent(in) :: mesh_dof
    real(dp),intent(in) ::  prq_solution(mesh_dof)
    integer,intent(in) :: depvar_at_elem(0:2,2,num_elems)
    integer,intent(in) :: depvar_at_node(num_nodes,0:2,2)
    !local variables
    integer :: np,ne,ny
    integer :: diagnostics_level


    character(len=60) :: sub_name
    sub_name = 'map_solution_to_mesh'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

      do  ne=1,num_elems
        ny=depvar_at_elem(1,1,ne)
        elem_field(ne_Qdot,ne)=prq_solution(ny)
      enddo !elems
      do np=1,num_nodes
        ny=depvar_at_node(np,0,1)
        node_field(nj_bv_press,np)=prq_solution(ny)
      enddo

      if(diagnostics_level.GT.1)then
        print *,"flow:"
        do ne=1,num_elems
          print *, "elem_field(",ne_Qdot,",",ne,")=",elem_field(ne_Qdot,ne)
        enddo
        print *,"pressure:"
        do np=1,num_nodes
          print *, "node_field(",nj_bv_press,",",np,")",node_field(nj_bv_press,np)
        enddo
      endif
    call enter_exit(sub_name,2)
end subroutine map_solution_to_mesh

!##############################################################################
!
subroutine map_flow_to_terminals
!*Description:* This subroutine maps the solution array to appropriate nodal and element fields
    use indices
    use arrays,only: elem_field,node_field,num_units,units,unit_field,elem_nodes
    use diagnostics, only: enter_exit,get_diagnostics_level
    integer :: nu,ne,np
    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'map_flow_to_terminals'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

	!num_units - number of terminal elements
    do nu=1,num_units
      ne=units(nu)!get terminal element
      np=elem_nodes(2,ne)
      unit_field(nu_perf,nu)=elem_field(ne_Qdot,ne)
      unit_field(nu_blood_press,nu)=node_field(nj_bv_press,np)
    enddo

    if(diagnostics_level.GT.1)then
      print *,"flow for terminal units:"
      do nu=1,num_units
        print *,"unit_field(",nu_perf,",",nu,")=",unit_field(nu_perf,nu)
      enddo
      print *,"pressure for terminal units:"
      do nu=1,num_units
        print *,"unit_field(",nu_blood_press,",",nu,")=",unit_field(nu_blood_press,nu)
      enddo
    endif

    call enter_exit(sub_name,2)
end subroutine map_flow_to_terminals
!
!##################################################################
!
subroutine get_variable_offset(depvar,mesh_dof,FIX,offset)
!*Description*: This subroutine returns the number of fixed variables in the depvar array that came before the input depvar
	integer, intent(in) :: depvar, mesh_dof
	logical, intent(in) :: FIX(mesh_dof)
    integer, intent(inout) :: offset

	integer :: depvar1

	offset = 0
	do depvar1 = 1,depvar
		if(FIX(depvar1))then
			offset = offset + 1
		endif
	enddo


end subroutine get_variable_offset
!
!##################################################################
!
subroutine calc_sparse_1dtree(bc_type,FIX,mesh_dof,depvar_at_elem,depvar_at_node,NonZeros,MatrixSize,&
SparseCol,SparseRow,SparseVal,RHS,prq_solution,update_resistance_entries,update_flow_nzz_row)
!*Description:* This subroutine sets up the system of equations to solve for pressure and flow. It populates the SparseCol,SparseRow, SparceVal and RHS vectors.
    use indices
    use arrays,only: dp,num_elems,elem_nodes,num_nodes,elems_at_node,elem_cnct,elem_field,num_arterial_elems,&
    umbilical_inlets,umbilical_outlets
    use diagnostics, only: enter_exit,get_diagnostics_level

	character(len=60),intent(in) :: bc_type
 	logical, intent(in) :: FIX(mesh_dof)
    integer, intent(in) :: mesh_dof
    integer,intent(in) :: depvar_at_elem(0:2,2,num_elems)
    integer,intent(in) :: depvar_at_node(num_nodes,0:2,2)
    integer, intent(in) :: NonZeros,MatrixSize

    integer, intent(inout) :: SparseCol(NonZeros)
    integer, intent(inout) :: SparseRow(MatrixSize+1)
    real(dp), intent(inout) :: SparseVal(NonZeros)
    real(dp), intent(inout) :: RHS(MatrixSize) !Right Hand Side of the linear system
    real(dp), intent(inout) :: prq_solution(mesh_dof)
    integer, intent(inout) :: update_resistance_entries(num_elems)
    integer, intent(inout) :: update_flow_nzz_row

	!local variables
    integer :: ne,nn,np,np1,np2,depvar,depvar1,depvar2,depvar3,flow_var,fixed_var_index,offset,nzz,&
      nzz_row,ne2,noelem2,ne3
    logical :: FlowBalancedNodes(num_nodes)
    logical :: NodePressureDone(num_nodes)
    logical :: ElementPressureEquationDone(num_elems)
    logical :: elem_found,one_node_balanced
    real(dp) :: flow_term
    integer :: inlet_counter
    logical :: unfixed_inlet = .False.
    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'calc_sparse_1dtree'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)


	inlet_counter = 0
	!Initialise matrices and indices
    SparseCol=0
    SparseRow=1
    SparseVal=0.0_dp
    RHS=0.0_dp

    nzz=1 !position in SparseCol and SparseVal
    nzz_row=1 !position in SparseRow
    FlowBalancedNodes = .FALSE. !.TRUE. for nodes which have had a conservation of flow equation done
    NodePressureDone = .FALSE.  !.TRUE. for nodes which have been processed
    ElementPressureEquationDone = .FALSE.
 	offset=0!variable position offset

    do ne=1,num_elems
  	  !look at pressure variables at each node
  	  do nn=1,2 !2 nodes in 1D element
		np=elem_nodes(nn,ne)
		depvar = depvar_at_node(np,1,1)
        if((.NOT.NodePressureDone(np)).AND.(.NOT.FIX(depvar)))then !check if this node is not fixed and hasn't already been processed (as nodes are shared between elements)
   	  		ne2=0
   	  		if(nn.EQ.1)then !first node of the element
   	  			ne2=ne! use the current element
   	  		elseif(nn.EQ.2)then !second node of the element
   	  			if((bc_type.EQ.'pressure').OR.(.NOT.ElementPressureEquationDone(ne)))then !if bc_type is pressure or element pressure equation for the current element hasn't been used
   	  				ne2=ne! use the current element
   	  			else
   	  				!look for another element connected to this node with pressure equation that hasn't been used
   	  				if (elems_at_node(np,0).GT.1)then
   	  					elem_found=.FALSE.
   	  					noelem2 = 1
   	  					do while ((.NOT.elem_found).AND.(noelem2.LE.elems_at_node(np,0)))
                			ne3=elems_at_node(np,noelem2)
                			if((ne3.NE.ne).AND.(.NOT.ElementPressureEquationDone(ne3)))then
                				ne2 = ne3
                				elem_found=.TRUE.
                			endif
                			noelem2 = noelem2 + 1
   	  					end do
   	  				endif
   	  			endif
   	  		endif
   	  		if((ne2.GT.0).and.(.NOT.ElementPressureEquationDone(ne2)))then
   	  			!do the pressure equation for element ne2
 				!pressure for node 1 - pressure for node 2 - resistance * flow at element ne2 = 0
				np1=elem_nodes(1,ne2)
	  			depvar1=depvar_at_node(np1,1,1) !pressure variable for first node
	  			np2=elem_nodes(2,ne2) !second node
      			depvar2=depvar_at_node(np2,1,1) !pressure variable for second node
	  			depvar3=depvar_at_elem(0,1,ne2) !flow variable for element
				if(FIX(depvar1))then !checking if pressure at 1st node is fixed
					!store known variable - inlet pressure
					RHS(nzz_row) = -prq_solution(depvar1)
				else
					!unknown variable -pressure for node 1
					call get_variable_offset(depvar1,mesh_dof,FIX,offset)
					SparseCol(nzz) = depvar1 - offset !variable number
					SparseVal(nzz)=1.0_dp !variable coefficient
					nzz=nzz+1 !next column
				endif
				if(FIX(depvar2))then !checking if pressure at 2nd node is fixed
	       			!store known variable - outlet pressure
					RHS(nzz_row) = prq_solution(depvar2)
				else
					!unknown variable - pressure for node 2
					call get_variable_offset(depvar2,mesh_dof,FIX,offset)
					SparseCol(nzz) = depvar2 - offset !variable number
					SparseVal(nzz)=-1.0_dp !variable coefficient
					nzz=nzz+1 !next column
				endif
				if(FIX(depvar3))then !checking if flow at element ne2 is fixed
					!store known variable - inlet flow * resistance for element	ne
					RHS(nzz_row) = prq_solution(depvar3)*elem_field(ne_resist,ne2)
					update_flow_nzz_row = nzz_row
				else
					!unknown flow
					call get_variable_offset(depvar3,mesh_dof,FIX,offset)
					SparseCol(nzz) = depvar3-offset !variable position in the unknown variable vector
					SparseVal(nzz)=-elem_field(ne_resist,ne2) !variable coefficient = resistance for element ne2
					update_resistance_entries(ne2) = nzz
					nzz=nzz+1 !next column
				endif
				nzz_row=nzz_row+1 !store next row position
	    			SparseRow(nzz_row)=nzz
   	  			NodePressureDone(np) = .TRUE.
   	  			ElementPressureEquationDone(ne2) = .TRUE.
   	  		endif
 	  	endif
  	  enddo !nn

  	  !look at flow variable for the element

	  flow_var = depvar_at_elem(0,1,ne)

	  if(.NOT.FIX(flow_var))then !don't do anything if flow is fixed
	  	  one_node_balanced = .FALSE.
	  	  !check if node 1 or node 2 are unbalanced
      	  do nn=1,2 !do flow balance for each element node
      	    np = elem_nodes(nn,ne)
            if((elems_at_node(np,0).GT.1).AND.(.NOT.FlowBalancedNodes(np)))then !if there is more than one element at a node and the node is not already flow balanced
          	  if((bc_type.EQ.'pressure').OR.((bc_type.EQ.'flow').AND.(.NOT.one_node_balanced)))then !do just one flow balance equation for bc_type flow
          		!go through each element connected to node np and add the conservation of flow equation for the elements
          		  do noelem2=1,elems_at_node(np,0)
              		    ne2=elems_at_node(np,noelem2)
              		    depvar=depvar_at_elem(1,1,ne2)
              		    flow_term = 0
              		    if(np.EQ.elem_nodes(2,ne2))then !end node
              			  flow_term = 1.0_dp
              		    elseif(np.EQ.elem_nodes(1,ne2))then !start node
              			  flow_term = -1.0_dp
              		    endif
              		    if(FIX(depvar))then
              			  RHS(nzz_row)=-prq_solution(depvar)*flow_term
              		    else
                		  !populate SparseCol and SparseVal
			  	          call get_variable_offset(depvar,mesh_dof,FIX,offset)
			  	          SparseCol(nzz) = depvar - offset
              			  SparseVal(nzz) = flow_term
              			  nzz = nzz + 1
              		    endif
			    enddo
				FlowBalancedNodes(np) = .TRUE.
				nzz_row=nzz_row+1 !store next row position
	    			SparseRow(nzz_row)=nzz
	    			one_node_balanced = .TRUE.
	    		endif !checking bc_type

 		    endif !flow at node np is unbalanced

 		  enddo !nn

	  	  !if flow balancing hasn't been done for any node for element ne and pressure equation hasn't already been done, do the pressure equation for the element
	  	  if((.NOT.one_node_balanced).AND.(.NOT.ElementPressureEquationDone(ne)))then

  	  		!do the pressure equation for element ne
 			!pressure for node 1 - pressure for node 2 - resistance * flow at element ne = 0
			np1=elem_nodes(1,ne)
	  		depvar1=depvar_at_node(np1,1,1) !pressure variable for first node
	  		np2=elem_nodes(2,ne) !second node
      		depvar2=depvar_at_node(np2,1,1) !pressure variable for second node

			!unknown variable -pressure for node 1
			call get_variable_offset(depvar1,mesh_dof,FIX,offset)
			SparseCol(nzz) = depvar1 - offset !variable number
			SparseVal(nzz)=1.0_dp !variable coefficient
			nzz=nzz+1 !next column

			if(FIX(depvar2))then !checking if pressure at 2nd node is fixed
	       		!store known variable - outlet pressure
				RHS(nzz_row) = prq_solution(depvar2)
			else
				!unknown variable - pressure for node 2
				call get_variable_offset(depvar2,mesh_dof,FIX,offset)
				SparseCol(nzz) = depvar2 - offset !variable number
				SparseVal(nzz)=-1.0_dp !variable coefficient
				nzz=nzz+1 !next column
			endif

			!unknown flow
			call get_variable_offset(flow_var,mesh_dof,FIX,offset)
			SparseCol(nzz) = flow_var-offset !variable position in the unknown variable vector
			SparseVal(nzz)=-elem_field(ne_resist,ne) !variable coefficient = resistance for element ne
			update_resistance_entries(ne) = nzz
			nzz=nzz+1 !next column

			nzz_row=nzz_row+1 !store next row position
	    		SparseRow(nzz_row)=nzz
   	  		ElementPressureEquationDone(ne) = .TRUE.
	  	  endif
	    endif
	    unfixed_inlet = .False.
	enddo !ne

    print *,"MatrixSize=",MatrixSize,nzz_row
    print *,"NonZeros=",NonZeros,nzz
    if(diagnostics_level.GT.1)then
        do nzz=1,NonZeros
  			print *,"SparseCol(",nzz,")",SparseCol(nzz)
  		enddo
  		do nzz_row=1,MatrixSize+1
   			print *,"SparseRow(",nzz_row,")",SparseRow(nzz_row)
  		enddo
  		do nzz=1,NonZeros
   			print *,"SparseVal(",nzz,")",SparseVal(nzz)
  		enddo
  		do nzz_row=1,MatrixSize
   			print *,"RHS(",nzz_row,")",RHS(nzz_row)
  		enddo
    endif
    call enter_exit(sub_name,2)
  end subroutine calc_sparse_1dtree
!##################################################
!

subroutine capillary_resistance(nelem,vessel_type,rheology_type,press_in,press_out,resistance,export_terminals)
    use indices
    use arrays
    use diagnostics
    use solve
    use other_consts
    integer, intent(in) :: nelem
    character(len=60), intent(in) :: rheology_type,vessel_type
    real(dp), intent(in) :: press_in
    real(dp), intent(in) :: press_out !When considering rigid vessels these can be dummy arguments, unless exporting, as resistance does not depend on pressure
    real(dp), intent(out) :: resistance
    logical, intent(in) :: export_terminals
    !Local variables to be passed in
    integer :: numconvolutes,numgens,numparallel
    !Actual local variables
    integer :: p_unknowns !Number of unknown pressures
    integer :: q_unknowns ! Number of unknown flows
    integer :: MatrixSize
    integer :: NonZeros
    integer, allocatable :: SparseCol(:)
    integer, allocatable ::SparseRow(:)
    real(dp),allocatable :: RHS(:)
    real(dp),allocatable :: Solution(:)
    real(dp), allocatable :: SolutionNew(:)
    real(dp),allocatable :: SparseVal(:)
    real(dp),allocatable :: original_resist(:)
    integer, allocatable :: update_resist(:,:)
    integer :: AllocateStatus
    integer :: SOLVER_FLAG
    integer :: ng,nc,i,i2,nnz,j,curgen,checkgen,nart,nv
    real(dp) :: int_length,int_radius,cap_length,cap_rad,seg_length,mu,R_seg,R_cap,divider
    real(dp) :: int_rad_ain,int_rad_aout,int_rad_vin, int_rad_vout,int_radius_gen
    integer :: update_mu,update_rad,count_its
    integer :: numseries,numparallel_cap
    real(dp) :: visc_factor,p_in,p_out,r_in,r_out,r_ave,h,elastance,err
    real(dp) :: total_art_surface_area,total_art_volume
    real(dp) :: total_vein_surface_area, total_vein_volume

    logical :: update_resistance = .False.
    logical :: converged = .False.
    character(len=60) :: sub_name
    integer :: diagnostics_level
    sub_name = 'capillary_resistance'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)
    if(capillary_model_type.eq.1)then
       mu=0.33600e-02_dp !; %viscosity
       resistance = 8.0_dp*mu*elem_field(ne_viscfact,nelem)*elem_field(ne_length,nelem)/ &
            (PI*((elem_field(ne_radius_in,nelem)+elem_field(ne_radius_out,nelem))/2.0_dp)**4.0_dp) !laminar resistance
    else
      numparallel = 1 !Number of convolute units in parallel !dummy variable of 1 but could be updated if needed
      !Actual variables
      numseries = 3 !Number of terminal villi in a row from a single mature intermediate villous

      numparallel_cap = num_parallel !Number of parallel capillaries in an imaged convolute (leiser)
      numconvolutes = num_convolutes!number of  terminal conduits in a single feeding vessel
      numgens = num_generations

      total_cap_volume =0.0_dp
      total_art_volume = 0.0_dp
      total_vein_volume = 0.0_dp
      total_cap_surface_area = 0.0_dp
      total_art_surface_area = 0.0_dp
      total_vein_surface_area = 0.0_dp

      int_length=1.5_dp!mm %Length of each intermediate villous
      int_radius=0.030_dp/2.0_dp!%radius of each intermediate villous
      seg_length=int_length/dble(numconvolutes) !; %lengh of each intermediate villous segment

      cap_length=3.0_dp/dble(numparallel_cap)!%mm %length of individual capillary
      cap_rad=(0.0144_dp)/2.0_dp ! %radius of individual capillary

      visc_factor = 1.0_dp !Initialise viscosity factor
      if(rheology_type.eq.'constant_visc')then
        mu=0.33600e-02_dp !; %viscosity
        update_mu = 0 !Do not need to update viscosity
      else if(rheology_type.eq.'pries_network')then
         mu=0.400e-02_dp !; %viscosity
         update_mu = 2 ! Need to update viscosity (iteratively)
      else if(rheology_type.eq.'pries_vessel')then
        mu=0.400e-02_dp !; %viscosity
        update_mu = 1 !Need to update viscosity, but not iteratively
      else !should have been caught earlier, but default to simplest case
        mu=0.33600e-02_dp !; %viscosity
        update_mu = 0 !Do not need to update viscosity
      endif
      if(vessel_type.eq.'rigid')then
        update_rad = 0
      else
        update_rad = 2
      endif

      nart = elem_cnct(-1,1,nelem) !capillary unit is downstream of a terminal unit
      nv =  elem_cnct(1,1,nelem) !vein is downstream of the capillary
      int_rad_ain= elem_field(ne_radius,nart) !mm Unstrained radius of inlet villous
      if(nv.gt.0)then
        int_rad_vin = elem_field(ne_radius,nv) !mm radius of ouutlet intermediate villous
      else
        int_rad_vin = elem_field(ne_radius,nart)*2.0_dp !mm radius of ouutlet intermediate villous
      end if
      int_rad_aout =  0.03_dp/2.0_dp ! mm radius of mature intermediate villous artery
      int_rad_vout = 0.03

      if(update_mu.eq.1) then
        call viscosity_from_radius(int_radius*1000.0_dp,0.45_dp,visc_factor)
      endif
      R_seg=(8.0_dp*mu*visc_factor*seg_length)/(pi*int_radius**4.0_dp)! %resistance of each intermediate villous segment
      if(update_mu.eq.1) then
        call viscosity_from_radius(cap_rad*1000.0_dp,0.45_dp,visc_factor)
      endif

      !In the Interface 2015 model the resistance of a "terminal capillary convolute) is defined anatomically by Leiser
      !Or, we can use the measured resistance of a "terminal capillary convolute" from Erlich
      if(capillary_model_type.eq.3)then
         R_cap =  5.6e5_dp*dble(numseries)/dble(numparallel)
         cap_resistance = 5.6e5_dp
      else if(capillary_model_type.eq.2)then
    	R_cap=(8.0_dp*mu*visc_factor*cap_length)/(pi*cap_rad**4.0_dp)*dble(numseries)/&
                (dble(numparallel_cap)*dble(numparallel))!%resistance of each capillary convolute segment
        cap_resistance = ((8.0_dp*mu*visc_factor*cap_length)/(pi*cap_rad**4.0_dp))/dble(numparallel_cap)
      endif

      p_unknowns = 2*numconvolutes*numgens
      q_unknowns = numgens*numconvolutes+1 !Only need to calculate convolute and total flows
      MatrixSize = p_unknowns+q_unknowns

    NonZeros =  q_unknowns*4 +  p_unknowns*q_unknowns/2 +p_unknowns*2-5
    !bttom row + capillaries, +  pyramid entries + bcs +pbalance

    allocate (SparseCol(NonZeros), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (SparseVal(NonZeros), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (SparseRow(MatrixSize+1), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (Solution(MatrixSize), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (RHS(MatrixSize), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"

    if((update_mu.eq.2).or.(update_rad.eq.2))then
       !Need to allocate memory for an iterative solution
       allocate (update_resist(NonZeros,3), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
        update_resistance = .True.
        update_resist = 0
       allocate (original_resist(NonZeros), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
        original_resist = 0.0_dp
        allocate (SolutionNew(MatrixSize), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    endif

    !Initialise solution matrix arrays
    SparseCol=0
    SparseRow=0
    SparseVal=0.0_dp
    Solution=0.0_dp
    RHS = 0.0_dp
    nnz=1
    !Put matrix together
    do ng = 1, numgens
      !Radius of branching villi at theis generation
      int_radius_gen = int_rad_ain + (int_rad_aout-int_rad_ain)/dble(numgens)*(dble(ng)-1)
      if(update_mu.eq.1) then
       call viscosity_from_radius(int_radius_gen*1000.0_dp,0.45_dp,visc_factor)
      endif
      !Resistance of a segment of branching villi between terminal convolute segments
      R_seg=(8.0_dp*mu*visc_factor*seg_length)/(pi*int_radius_gen**4.0_dp)! %resistance of each intermediate villous segment
      !Add volume and surface area to this segment
      total_art_volume = total_art_volume +  (PI*int_radius_gen**2.0_dp*seg_length*dble(numconvolutes))*(2.0_dp**dble(ng))
      total_art_surface_area = total_art_surface_area + PI*2.0_dp*int_radius_gen*seg_length*dble(numconvolutes)&
      	*(2.0_dp**dble(ng))
      !Add corresponding "capillary" or terminal conduit volume and surface area to this segment
      if(capillary_model_type.eq.3)then
        total_cap_volume = total_cap_volume +0.65e-3_dp*dble(numconvolutes)*dble(numseries)*dble(numparallel)&
                *(2.0_dp**dble(ng))
        total_cap_surface_area = total_cap_surface_area + 0.135_dp*dble(numconvolutes)*dble(numseries)* &
              dble(numparallel)*(2.0_dp**dble(ng))
      else if(capillary_model_type.eq.2)then
        total_cap_volume = total_cap_volume +   PI*cap_rad**2.0_dp*dble(cap_length)*dble(numconvolutes)*&
                dble(numseries)*dble(numparallel)*dble(numparallel_cap)*(2.0_dp**(ng))
        total_cap_surface_area = 2.0_dp*PI*cap_rad*cap_length*dble(numconvolutes)*dble(numseries)*dble(numparallel)* &
                dble(numparallel_cap)*(2.0_dp**dble(ng))
      endif
      !Add matrix entries for each convolute within this generation
      do nc = 1,numconvolutes
        i = (ng-1)*numconvolutes + nc !row number
        !outward branches (arteries)
        if((nc.eq.1).and.(ng.eq.1))then !Inlet -pressure out -QR =-Pressure in
          RHS(i) = -press_in
          SparseRow(i) = 1
          SparseCol(nnz) = i
          SparseVal(nnz) = -1.0_dp
          nnz=nnz+1
          SparseCol(nnz) = MatrixSize
          SparseVal(nnz) = -R_seg/2.0_dp !divide by 2 because bifurcating here
          if(update_resistance)then
            update_resist(nnz,1) = 1 !indicating this one needs to be updated and is an artery
            update_resist(nnz,2) = 1 !power for resistance division
            update_resist(nnz,3) = i !segment that this resistance is associated with
            original_resist(nnz) = SparseVal(nnz)
          endif
          nnz = nnz+1
        else
          SparseCol(nnz) = i-1
          SparseVal(nnz) = 1.0_dp
          nnz = nnz+1
          SparseCol(nnz) = i
          SparseVal(nnz) = -1.0_dp
          nnz=nnz+1
          divider = 2.0_dp**(ng-1) !FLOW DIVIDES BY 2
          curgen = ng
          checkgen = 1
          do j = 1,i-1
            SparseCol(nnz) = p_unknowns + j
            SparseVal(nnz) = R_seg/divider
            if(update_resistance)then
              update_resist(nnz,1) = 1 !indicating this one needs to be updated and is an artery
              update_resist(nnz,2) = ng !current generation
              update_resist(nnz,3) = j !segment that this resistance is associated with
              original_resist(nnz) = SparseVal(nnz)
            endif
            nnz = nnz+1
            if(j.eq.(checkgen)*numconvolutes)then
              checkgen = checkgen + 1
              curgen = curgen - 1
              divider = 2.0_dp**(curgen - 1)
            endif
          enddo
          SparseCol(nnz) = MatrixSize
          SparseVal(nnz) = -R_seg/(2.0_dp**ng)
          if(update_resistance)then
            update_resist(nnz,1) = 1 !indicating this one needs to be updated and is an artery
            update_resist(nnz,2) = ng !current generation
            update_resist(nnz,3) = i !segment that this resistance is associated with
            original_resist(nnz) = SparseVal(nnz)
          endif
          nnz = nnz+1
        endif
        SparseRow(i+1) = nnz
      enddo
    enddo
    !Repeat for veins
    do ng = 1, numgens

      int_radius_gen = int_rad_vin + (int_rad_vout-int_rad_vin)/dble(numgens)*(dble(ng)-1.0_dp)

      if(update_mu.eq.1) then
        call viscosity_from_radius(int_radius_gen*1000.0_dp,0.45_dp,visc_factor)
      endif
      R_seg=(8.0_dp*mu*visc_factor*seg_length)/(pi*int_radius_gen**4.0_dp)! %resistance of each intermediate villous segment

      total_vein_volume = total_vein_volume +  PI*int_radius_gen**2.0_dp*seg_length*dble(numconvolutes)*&
              (2.0_dp**dble(ng))
      total_vein_surface_area = total_vein_surface_area + PI*2.0_dp*int_radius_gen*seg_length*dble(numconvolutes)&
      	*(2.0_dp**dble(ng))
      do nc = 1,numconvolutes
        i = (ng-1)*numconvolutes + nc !pointer to arteries row number
        !inward branches (veins)
        i2 =p_unknowns/2 + (ng-1)*numconvolutes + nc
        SparseCol(nnz) = i2
        SparseVal(nnz) = 1.0_dp
        nnz = nnz+1
        if(i2.eq.p_unknowns)then
          RHS(i2) = press_out
          SparseCol(nnz) = MatrixSize
          SparseVal(nnz) = -R_seg/2.0_dp
          if(update_resistance)then
            update_resist(nnz,1) = 2 !indicating this one needs to be updated and is a vein
            update_resist(nnz,2) = numgens - ng + 1 !current generation
            update_resist(nnz,3) = i2 !segment that this resistance is associated with
            original_resist(nnz) = SparseVal(nnz)
          endif
          nnz = nnz + 1
        else
          SparseCol(nnz) = i2+1
          SparseVal(nnz) = -1.0_dp
          nnz = nnz+1
          checkgen = 1
          curgen = ceiling((dble(p_unknowns) - dble(i2))/dble(numconvolutes))
          divider = 2.0_dp**(curgen-1)
          do j = 1,(p_unknowns - i2)
            if(j.eq.(p_unknowns - i2))then
              SparseCol(nnz) = p_unknowns+j
              SparseVal(nnz) = R_seg
              if(update_resistance)then
                update_resist(nnz,1) = 2 !indicating this one needs to be updated and is a vein
                update_resist(nnz,2) = numgens - ng + 1 !current generation
                update_resist(nnz,3) = i2 !segment that this resistance is associated with
                original_resist(nnz) = SparseVal(nnz)
              endif
              nnz = nnz+1
            else
              SparseCol(nnz) = p_unknowns + j
              SparseVal(nnz) = R_seg/divider
              if(update_resistance)then
                update_resist(nnz,1) = 2 !indicating this one needs to be updated and is a vein
                update_resist(nnz,2) = numgens - ng + 1 !current generation
                update_resist(nnz,3) = i2 !segment that this resistance is associated with
                original_resist(nnz) = SparseVal(nnz)
              endif
              nnz = nnz + 1
            endif
            if(j.eq.checkgen*numconvolutes)then
              checkgen = checkgen + 1
              curgen = curgen - 1
              divider = 2.0_dp**(curgen-1.0_dp)
            endif
          enddo
          curgen = ceiling((dble(p_unknowns) - dble(i2))/dble(numconvolutes))
          SparseCol(nnz) = MatrixSize
          SparseVal(nnz) = -R_seg/(2.0_dp**curgen)
          if(update_resistance)then
            update_resist(nnz,1) = 2 !indicating this one needs to be updated and is a vein
            update_resist(nnz,2) = numgens - ng + 1!current generation
            update_resist(nnz,3) = i2 !segment that this resistance is associated with
            original_resist(nnz) = SparseVal(nnz)
          endif
          nnz = nnz+1
        endif
        SparseRow(i2+1) = nnz
      enddo
     enddo
    do ng = 1, numgens
      do nc = 1,numconvolutes
        !Capillary segments
        i = (ng-1)*numconvolutes + nc !pointer to arteries row number
        SparseCol(nnz) = i
        SparseVal(nnz) = 1.0_dp
        nnz = nnz+1
        SparseCol(nnz) = p_unknowns + 1 - i
        SparseVal(nnz) = -1.0_dp
        nnz = nnz+1
        SparseCol(nnz) = p_unknowns+i
        SparseVal(nnz) = -R_cap!/2.0_dp
        nnz = nnz+1
        SparseRow(p_unknowns+i+1) = nnz
      enddo
    enddo

    do ng = 1, numgens
      do nc = 1,numconvolutes
        !Conservation of flow
        i = (ng-1)*numconvolutes + nc !pointer to arteries row number
        SparseCol(nnz) = p_unknowns+i
        SparseVal(nnz) = 2.0_dp**ng
        if(update_resistance)then
            update_resist(nnz,1) = 3 !indicating this one needs to be updated and is a vein
            update_resist(nnz,2) = ng!current generation
            update_resist(nnz,3) = i !segment that this resistance is associated with
            original_resist(nnz) = SparseVal(nnz)
        endif
        nnz = nnz+1
      enddo
    enddo
    SparseCol(nnz) = MatrixSize
    SparseVal(nnz) = -1.0_dp
    nnz = nnz+1
    SparseRow(MatrixSize+1) = nnz

    call pmgmres_ilu_cr(MatrixSize, NonZeros, SparseRow, SparseCol, SparseVal, Solution, &
           RHS, 500, 500,1.d-5,1.d-4, SOLVER_FLAG)

    !!!!!!!!!!!!!!!ONLY IF NOT RIGID!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !We now have pressures - which we can use to update radii and resistances and re-solve
    h=0.8_dp
    elastance = 5.0e5_dp!
    if(update_resistance)then
       err = 0.0_dp
       SolutionNew = Solution
       count_its = 0
       converged = .False.
      do while (.not.converged)
       do nnz = 1,NonZeros
         if (update_resist(nnz,1).eq.1)then!artery
           if(update_resist(nnz,3).eq.1)then
             p_in = press_in -11.0_dp*133.0_dp
           else
             p_in  = Solution(update_resist(nnz,3)-1)-11.0_dp*133.0_dp
           endif
           p_out = Solution(update_resist(nnz,3))-11.0_dp*133.0_dp
           int_radius_gen = int_rad_ain + (int_rad_aout-int_rad_ain)/dble(numgens)*(dble(update_resist(nnz,2))-1.0_dp)

           r_in = int_radius_gen+3.0_dp*int_radius_gen**2*p_in/(4.0_dp*elastance*h)
           r_out = int_radius_gen+3.0_dp*int_radius_gen**2*p_out/(4.0_dp*elastance*h)
           r_ave = (r_in + r_out)/2.0_dp
           SparseVal(nnz) = original_resist(nnz)*(int_radius_gen**4.0_dp)/(r_ave**4.0_dp)
         elseif (update_resist(nnz,1).eq.2)then!vein
           if(update_resist(nnz,3).eq.p_unknowns)then
             p_out = press_out -11.0_dp*133.0_dp
           else
             p_out  = Solution(update_resist(nnz,3)+1)-11.0_dp*133.0_dp
           endif
           p_in = Solution(update_resist(nnz,3))-11.0_dp*133.0_dp
           int_radius_gen = int_rad_vin + (int_rad_vout-int_rad_vin)/dble(numgens)*(dble(update_resist(nnz,2))-1.0_dp)

           r_in = int_radius_gen+3.0_dp*int_radius_gen**2*p_in/(4.0_dp*elastance*h)
           r_out = int_radius_gen+3.0_dp*int_radius_gen**2*p_out/(4.0_dp*elastance*h)
           r_ave = (r_in + r_out)/2.0_dp
           SparseVal(nnz) = original_resist(nnz)*(int_radius_gen**4.0_dp)/(r_ave**4.0_dp)
         elseif (update_resist(nnz,1).eq.3)then!capillary
           if(capillary_model_type.eq.2)then !We have a simple branching model and can update capillary resistamc
             p_in = Solution(update_resist(nnz,3))-11.0_dp*133.0_dp
             p_out = Solution(p_unknowns-update_resist(nnz,3)+1)-11.0_dp*133.0_dp
             r_in = cap_rad+3.0_dp*cap_rad**2*p_in/(4.0_dp*elastance*h)
             r_out = cap_rad+3.0_dp*cap_rad**2*p_out/(4.0_dp*elastance*h)
             r_ave = (r_in + r_out)/2.0_dp
             SparseVal(nnz) = original_resist(nnz)*(cap_rad**4.0_dp)/(r_ave**4.0_dp)
           else!We have an anatomically derived Erlich Model
               SparseVal(nnz) = original_resist(nnz)
           end if
         endif
        enddo
       ! Solve system
       call pmgmres_ilu_cr(MatrixSize, NonZeros,SparseRow, SparseCol, SparseVal, SolutionNew, &
               RHS, 500, 500,1.d-5,1.d-4, SOLVER_FLAG)
       do i = 1,MatrixSize
         err = err + (SolutionNew(i) - Solution(i))**2.0_dp/Solution(i)**2.0_dp
       enddo
       err = err/MatrixSize
       Solution = SolutionNew
       if(err.lt.1.0e-6_dp)then
         converged = .True.
       else if (count_its.gt.20)then
         converged = .True.
       endif
       count_its = count_its + 1
       enddo !while
    endif
    !!Total resistance of capillary system
    resistance = (press_in-press_out)/(Solution(MatrixSize))
    !!Effective length of capillary system
    terminal_length = resistance*(PI*(0.03_dp)**4.0_dp)/(8.0_dp*mu)

    if(export_terminals)then
       !write to file
        WRITE(10,&
              '(I6,X,2(F9.3,X),F9.5,X)', Advance = 'No' )&!,I6,X,4(F8.2,X),4(F8.5,X),&2(F10.2,X),3(F8.4,X),I6,X,2(F10.5,X),2(F8.4,X),&(F10.2,X))')&
         nelem,press_in,press_out,Solution(MatrixSize)
         do i = 1,MatrixSize-1
            if(i.lt.p_unknowns)then
              WRITE(10,'(F9.3,X)', Advance = 'No') Solution(i)
            else
              WRITE(10,'(F9.5,X)', Advance = 'No') Solution(i)
            endif
         enddo
         write(10, '(2(I6,X))') numgens,numconvolutes
    endif

    deallocate(SparseCol,STAT=AllocateStatus)
    deallocate(SparseVal,STAT=AllocateStatus)
    deallocate(SparseRow,STAT=AllocateStatus)
    deallocate(Solution,STAT=AllocateStatus)
    deallocate(RHS,STAT=AllocateStatus)
    if(update_resistance)then
       deallocate (update_resist, STAT = AllocateStatus)
       deallocate (original_resist, STAT = AllocateStatus)
       deallocate (SolutionNew, STAT = AllocateStatus)
    endif
    !Update fields for volume and surface area (multiply by two due to branch on entry into capillary system)
    elem_field(ne_vol,nelem) = total_cap_volume
    elem_field(ne_sa,nelem) = total_cap_surface_area
    elem_field(ne_artvol,nelem) = total_art_volume
    elem_field(ne_artsa,nelem) = total_art_surface_area
    elem_field(ne_veinvol,nelem) = total_vein_volume
    elem_field(ne_veinsa,nelem) = total_vein_surface_area
    elem_field(ne_length,nelem) = terminal_length

    end if

    call enter_exit(sub_name,2)
end subroutine capillary_resistance

!
!##################################################################
!

!
!#################################################################3
!
subroutine viscosity_from_radius(radius,hb,visc_factor)
    use indices
    use arrays
    use diagnostics

    real(dp), intent(in) :: radius
    real(dp), intent(in) :: hb
    real(dp), intent(inout) :: visc_factor


    real(dp) :: controlhb
    real(dp) :: beta

    character(len=60) :: sub_name
    integer :: diagnostics_level

  sub_name = 'map_flow_to_terminals'
  call enter_exit(sub_name,1)
  call get_diagnostics_level(diagnostics_level)


  controlhb = 0.45_dp
  beta = 4.0_dp/(1+exp(-0.0593_dp*(2.0_dp*radius-6.74_dp)))
  visc_factor = (1.0_dp + (exp(hb*beta)-1.0_dp)/(exp(controlhb*beta)-1.0_dp) &
          *(110.0_dp*exp(-2.848_dp*radius) + 3.0_dp - 3.45_dp*exp(-0.07_dp*radius)))/4.0_dp

  call enter_exit(sub_name,2)
end subroutine viscosity_from_radius



end module pressure_resistance_flow

