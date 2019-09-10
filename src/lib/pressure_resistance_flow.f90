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
  public evaluate_prq, calculate_stats
contains
!###################################################################################
!
subroutine evaluate_prq(mesh_type,bc_type,inlet_flow,inlet_pressure,outlet_pressure)
!*Description:* Solves for pressure and flow in a rigid or compliant tree structure  
! Model Types:                                                                     
! mesh_type: can be 'simple_tree' or 'full_plus_tube'. Simple_tree is the input arterial tree 
! without any special features at the terminal level
! 'full_plus_tube' creates a matching venous mesh and has arteries and veins connected by
! capillary units (capillaries are just tubes represented by an element)
!                                                                                  
! boundary condition type (bc_type): pressure or flow
    use indices
    use arrays,only: dp,num_elems,num_nodes,elem_field,elem_nodes,elem_cnct,node_xyz
    use diagnostics, only: enter_exit,get_diagnostics_level
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_PRQ" :: EVALUATE_PRQ
  
    character(len=60), intent(in) :: mesh_type,bc_type
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
    integer :: num_vars,NonZeros,MatrixSize
    integer :: AllocateStatus
    integer :: diagnostics_level

    real(dp), allocatable :: prq_solution(:),solver_solution(:)
    real(dp) :: viscosity,density,inletbc,outletbc,gamma,total_resistance
    logical, allocatable :: FIX(:)
    logical :: ADD=.FALSE.,CONVERGED=.FALSE.
    character(len=60) :: sub_name
    integer :: no,depvar,nz,ne,SOLVER_FLAG,ne0,ne1,nj,i

    sub_name = 'evaluate_prq'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

if((mesh_type.NE.'full_plus_tube').AND.(mesh_type.NE.'simple_tree'))then
	print *,"unsupported mesh_type",mesh_type
	call exit(1)
endif

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

if(diagnostics_level.GT.1)then
	print *, "mesh_type=",mesh_type
	print *, "bc_type=",bc_type
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

!!---------PHYSICAL PARAMETERS-----------
density=0.10500e-02_dp !kg/cm3 !density:fluid density
viscosity=0.33600e-02_dp !Pa.s !viscosity: fluid viscosity
gamma = 0.327_dp !=1.85/(4*sqrt(2)) !gamma:Pedley correction factor

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
        
!! Calculate sparsity structure for solution matrices
    !Determine size of and allocate solution vectors/matrices
    call calc_sparse_size(mesh_dof,FIX,depvar_at_elem,MatrixSize,NonZeros)

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
        prq_solution)    
 
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
          
!! ----CALL SOLVER----
   call pmgmres_ilu_cr(MatrixSize, NonZeros, SparseRow, SparseCol, SparseVal, &
         solver_solution, RHS, 500, 500,1.d-5,1.d-4,SOLVER_FLAG)
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
    
!need to write solution to element/nodal fields for export
    call map_solution_to_mesh(prq_solution,depvar_at_elem,depvar_at_node,mesh_dof)
    !NEED TO UPDATE TERMINAL SOLUTION HERE. LOOP THO' UNITS AND TAKE FLOW AND PRESSURE AT TERMINALS
    call map_flow_to_terminals

    deallocate (mesh_from_depvar, STAT = AllocateStatus)
    deallocate (depvar_at_elem, STAT = AllocateStatus)
    deallocate (depvar_at_node, STAT = AllocateStatus)
    deallocate (prq_solution, STAT = AllocateStatus)
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
    integer :: nonode,np,ne,ny1,nj
    character(len=60) :: sub_name

  sub_name = 'boundary_conditions'
  call enter_exit(sub_name,1)
  if(.NOT.ADD)THEN
     ! Initial values
     FIX(1:mesh_dof)=.FALSE.
     prq_solution = 0
     ! Fixed boundary conditions  
     ! These are inlet BCs, apply to all inlet BCs (there should only be one)
     do ne=1,num_elems
        !ne=elems(noelem)
        if (elem_cnct(-1,0,ne) == 0) THEN !Entry element
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
    real(dp) :: resistance,zeta
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
       ! element Poiseuille (laminar) resistance in units of Pa.s.mm-3        
       resistance = 8.d0*viscosity*elem_field(ne_length,ne)/ &
            (PI*elem_field(ne_radius,ne)**4) !laminar resistance
       elem_field(ne_resist,ne) = resistance
       if(diagnostics_level.GT.1)then
       		print *,"TESTING RESISTANCE: element",ne,"resistance",elem_field(ne_resist,ne),"radius",elem_field(ne_radius,ne)
       endif
    enddo 

    call enter_exit(sub_name,2)
  end subroutine calculate_resistance
!
!##################################################################
!
subroutine calculate_stats()
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
    use indices
    use arrays,only: num_arterial_elems,dp,num_elems, &
                  elem_field,num_units,units,elem_cnct,elem_ordrs, &
                  num_conv,num_conv_gen,cap_resistance,terminal_resistance, &
                  terminal_length,total_vasc_resistance, &
                  cap_radius
    use other_consts, only: PI,MAX_FILENAME_LEN
    use diagnostics, only: enter_exit,get_diagnostics_level
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CALCULATE_STATS" :: CALCULATE_STATS
  
  !local variables
    character(len=MAX_FILENAME_LEN) :: FLOW_GEN_FILE
    integer :: ne,nu,nc,order,max_strahler,no_branches, &
               ne_order,branch,ven_elems,max_gen
    real(dp) :: arterial_vasc_volume, total_vasc_volume, capillary_volume, &
                venous_vasc_volume, single_cap_surface_area, total_cap_surface_area, &
                mean_diameter,std_diameter,total_resistance,std_terminal_flow, &
                cof_var_terminal_flow,mean_terminal_flow,small_vessel_volume, diameter, &
                image_voxel_size,single_capillary_volume, cap_length,cof_var_cap_flow
    integer :: strahler_orders(num_elems)
    integer :: generations(num_elems)
    integer :: capillaries(num_units)
    real(dp),allocatable :: diameter_by_strahler(:,:)
    real(dp),allocatable :: art_diameter_by_strahler(:,:)
    real(dp),allocatable :: ven_diameter_by_strahler(:,:)
    real(dp),allocatable :: terminal_flow_by_gen(:,:)
    integer,allocatable :: branch_count(:)
    integer,allocatable :: gen_branch_count(:)
    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'calculate_stats'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)


    image_voxel_size = 0.1165 !mm

  !calculate arterial vascular volume
   arterial_vasc_volume = 0
   do ne=1,num_arterial_elems
      arterial_vasc_volume = arterial_vasc_volume + elem_field(ne_vol,ne)
   enddo
   print *, "Arterial elems count",num_arterial_elems
   print *, "Arterial vascular volume (cm**3) = ",arterial_vasc_volume/1000

   capillary_volume = 0
   cap_length=3_dp
   single_capillary_volume = PI * cap_radius**2 * cap_length 
   !single capillary volume * number of branching generations * number of terminal units
   capillary_volume = single_capillary_volume * num_conv_gen * num_units
 
   print *, "Capillary unit count =",num_units
   print *, "Capillary volume (cm**3) = ",capillary_volume/1000

   !calculate capillary unit surface area

   !surface area = 2*PI*radius*length
   single_cap_surface_area = 2 * PI * cap_radius * cap_length
   !single capillary surface area * number of branching generations * number of terminal units
   total_cap_surface_area = single_cap_surface_area * num_conv_gen * num_units     
   print *, "Total capillary surface area (cm**2) = ", total_cap_surface_area/100    

   total_vasc_volume = 0
   do ne = 1,num_elems
      total_vasc_volume = total_vasc_volume + elem_field(ne_vol,ne)
   enddo
   print *, "Total vascular volume (cm**3) = ",total_vasc_volume/1000

   venous_vasc_volume = total_vasc_volume - arterial_vasc_volume - capillary_volume  
   print *, "Venous vascular volume (cm**3) = ",venous_vasc_volume/1000

   !volume of vessels with diameter smaller than 0.1165 mm (resolution of images)
   small_vessel_volume = 0
   do ne = 1,num_elems
      diameter = elem_field(ne_radius,ne) * 2
      if(diameter.LT.image_voxel_size)then
         small_vessel_volume = small_vessel_volume + elem_field(ne_vol,ne)
      endif
   enddo
   print *, "Volume of vessels with diameter smaller than 0.1165 mm (cm**3) = ",small_vessel_volume/1000

   !Print capillary convolute number, generations and total length of capillaries
   print *, "Number of capillary convolutes per generation = ",num_conv
   print *, "Number of capillary generations per terminal unit = ",num_conv_gen
   print *, "Resistance of capillary conduits (Pa.s/mm**3)=",cap_resistance
   print *, "Resistance of all generations of capillaries per terminal unit (Pa.s/mm**3)=",terminal_resistance
   print *, "Effective length of each terminal unit (mm)",terminal_length

   !total vascular resistance
   print *, "Total vascular resistance (Pa.s/mm**3) = ",total_vasc_resistance

   !mean, min, max and std of branch diameter by Strahler order

   strahler_orders(:) = elem_ordrs(no_sord, :)
   max_strahler = maxval(strahler_orders)
   allocate(diameter_by_strahler(max_strahler,num_elems))
   allocate(branch_count(max_strahler))
   diameter_by_strahler = 0
   branch_count = 0

   do ne=1,num_elems
      ne_order = strahler_orders(ne)
      no_branches = branch_count(ne_order)
      no_branches = no_branches + 1
      diameter_by_strahler(ne_order,no_branches) = elem_field(ne_radius,ne) * 2
      branch_count(ne_order) = no_branches
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


   !Arterial vessel diameter by Strahler order
   allocate(art_diameter_by_strahler(max_strahler,num_arterial_elems))
   art_diameter_by_strahler = 0
   branch_count = 0
   do ne=1,num_arterial_elems
      ne_order = strahler_orders(ne)
      no_branches = branch_count(ne_order)
      no_branches = no_branches + 1
      art_diameter_by_strahler(ne_order,no_branches) = elem_field(ne_radius,ne) * 2
      branch_count(ne_order) = no_branches
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
   
   !populate capillaries array
   do nu=1,num_units
      !get the element below a terminal branch - this is the capillary unit
      ne =units(nu) !Get a terminal unit   
      nc = elem_cnct(1,1,ne) !capillary unit is downstream of a terminal unit
      capillaries(nu) = nc
   enddo

   !all elements - arterial elements - number of capillaries (the same as number of terminal units)
   ven_elems = num_elems-num_arterial_elems-num_units 
   if (ven_elems.GT.0)then
      allocate(ven_diameter_by_strahler(max_strahler,ven_elems)) 
      ven_diameter_by_strahler = 0
      branch_count = 0
      do ne=num_arterial_elems+1,num_elems
         if(ALL(capillaries.NE.ne))then !if the element is not a capillary
            ne_order = strahler_orders(ne)
            no_branches = branch_count(ne_order)
            no_branches = no_branches + 1
            ven_diameter_by_strahler(ne_order,no_branches) = elem_field(ne_radius,ne) * 2
            branch_count(ne_order) = no_branches
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
  

  !coefficient of variation for terminal flow
   !standard deviation of flow devided by the mean flow
   mean_terminal_flow = 0
   do nu=1,num_units
      ne = units(nu)
      mean_terminal_flow = mean_terminal_flow + elem_field(ne_Qdot,ne)
   enddo
   mean_terminal_flow = mean_terminal_flow/num_units
   std_terminal_flow = 0
   do nu=1,num_units
      ne = units(nu)
      std_terminal_flow = std_terminal_flow + (elem_field(ne_Qdot,ne) - mean_terminal_flow)**2
   enddo
   std_terminal_flow = std_terminal_flow/num_units
   std_terminal_flow = SQRT(std_terminal_flow)
   cof_var_terminal_flow = std_terminal_flow/mean_terminal_flow
   print *, "Coefficient of variation for terminal flow (%) = ", cof_var_terminal_flow * 100


  
   !terminal flow by generation
   generations(:) = elem_ordrs(no_gen, :)
   max_gen = maxval(generations)
   allocate(terminal_flow_by_gen(max_gen,num_units))   
   allocate(gen_branch_count(max_gen))
   terminal_flow_by_gen = 0
   gen_branch_count = 0

  !print all terminal flows and their corresponding generations to a file 
   FLOW_GEN_FILE = "Output/terminal flow per generation.csv"
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


   call enter_exit(sub_name,2)

end subroutine calculate_stats
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
!*Descripton:* This subroutine calculates the total resistance of a tree
    use indices
    use arrays,only: dp,num_elems,elem_cnct,elem_field,&
                     total_vasc_resistance
    use diagnostics, only: enter_exit
    character(len=60) :: sub_name
!local variables
    real(dp), intent(out) :: resistance
    real(dp) :: invres,elem_res(num_elems)
    integer :: num2,ne,ne2

    sub_name = 'tree_resistance'
    call enter_exit(sub_name,1)

    elem_res(1:num_elems)=elem_field(ne_resist,1:num_elems)
    do ne=num_elems,1,-1
       !ne=elems(num)
      invres=0.0_dp
      do num2=1,elem_cnct(1,0,ne)
         ne2=elem_cnct(1,num2,ne)
         invres=invres+1.0_dp/elem_res(ne2) !resistance in parallel, for daughter branches
      enddo
      if(elem_cnct(1,0,ne).gt.0)then 
        elem_res(ne)=elem_res(ne)+1.0_dp/invres !resistance in a series
       endif
    enddo
    resistance=elem_res(1)
    total_vasc_resistance = resistance

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
subroutine calc_sparse_size(mesh_dof,FIX,depvar_at_elem,MatrixSize,NonZeros)
!*Description:* This subroutine calculates sparsity sizes

    use diagnostics, only: enter_exit,get_diagnostics_level
    use arrays,only: num_elems,elem_nodes,num_nodes,elems_at_node
       
    integer, intent(in) :: mesh_dof
    logical, intent(in) :: FIX(mesh_dof)
    integer,intent(in) :: depvar_at_elem(0:2,2,num_elems)
    integer,intent(inout) :: MatrixSize
    integer,intent(inout) :: NonZeros
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
SparseCol,SparseRow,SparseVal,RHS,prq_solution)
!*Description:* This subroutine sets up the system of equations to solve for pressure and flow. It populates the SparseCol,SparseRow, SparceVal and RHS vectors.
    use indices
    use arrays,only: dp,num_elems,elem_nodes,num_nodes,elems_at_node,elem_cnct,elem_field
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

	!local variables
    integer :: ne,nn,np,np1,np2,depvar,depvar1,depvar2,depvar3,flow_var,fixed_var_index,offset,nzz,&
      nzz_row,ne2,noelem2,ne3
    logical :: FlowBalancedNodes(num_nodes)
    logical :: NodePressureDone(num_nodes)
    logical :: ElementPressureEquationDone(num_elems)
    logical :: elem_found,one_node_balanced
    real(dp) :: flow_term
    character(len=60) :: sub_name
    integer :: diagnostics_level
    
    sub_name = 'calc_sparse_1dtree'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)
	
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
   	  		if(ne2.GT.0)then 
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
				else
					!unknown flow
					call get_variable_offset(depvar3,mesh_dof,FIX,offset)			
					SparseCol(nzz) = depvar3-offset !variable position in the unknown variable vector
					SparseVal(nzz)=-elem_field(ne_resist,ne2) !variable coefficient = resistance for element ne2
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
			nzz=nzz+1 !next column	 				
					
			nzz_row=nzz_row+1 !store next row position
	    		SparseRow(nzz_row)=nzz 			
   	  		ElementPressureEquationDone(ne) = .TRUE.	  		  
	  	endif	  
	  endif
	enddo !ne
	if(diagnostics_level.GT.1)then    
    		print *,"MatrixSize=",MatrixSize
    		print *,"NonZeros=",NonZeros
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


end module pressure_resistance_flow

