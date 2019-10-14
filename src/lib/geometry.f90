module geometry
!*Description:* This module handles all geometry read/write/generation.
!
  use other_consts
  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private
  public add_matching_mesh
  public append_units
  public calc_capillary_unit_length
  public define_1d_elements
  public define_node_geometry
  public define_rad_from_file
  public define_rad_from_geom
  public define_ven_rad_from_art
  public element_connectivity_1d
  public evaluate_ordering
  public get_final_real
  public update_1d_elem_field

contains


!
!###################################################################################
!
  subroutine add_matching_mesh(umbilical_elem_option_in, umbilical_element_numbers)
  !*Description:* adds a matching venous mesh to an arterial mesh.
  ! Two options are availabe: 1) 'same_as_arterial' - the venous vessels
  ! are an exact copy of the arterial vessels. 2) 'single_umbilical_vein' -
  ! a file containing umbilical arterial elements in separate lines must be
  ! given. The umbilical arterial vessels will not be copied. 
  ! A single umbilical vein will be created instead.
 
    use arrays,only: dp,elems,elem_cnct,elem_direction,elem_field,&
         elem_nodes,elem_ordrs,elem_symmetry,elems_at_node,&
         nodes,node_xyz,num_elems,&
         num_nodes,num_units,units,num_arterial_elems,umbilical_inlets, &
         umbilical_outlets, art_ven_elem_map
    use indices
    use other_consts,only: PI,MAX_STRING_LEN
    use diagnostics, only: enter_exit,get_diagnostics_level
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ADD_MATCHING_MESH" :: ADD_MATCHING_MESH 
    character(len=MAX_STRING_LEN),intent(in) ::  umbilical_elem_option_in  !argument controlling how the umbilical 
    !venous elements are created, 'same_as_arterial' or 'single_umbilical_vein'
    integer, intent(in), optional :: umbilical_element_numbers(:)  !list containing umbilical artery element numbers
    !Parameters to become inputs
    real(dp) :: offset(3)
    logical :: REVERSE=.TRUE.
    !local variables
    character(LEN=MAX_STRING_LEN) :: umbilical_elem_option
    integer :: num_nodes_new,num_elems_new,ne,ne_global,np,np_global,np0,nonode,np_m, &
            downstream_elem, node_counter, elem_counter, indx, elem_at_node, indx2, &
            nj,ne_m,noelem,ne0,n,nindex,ne1,noelem0,nu,cap_conns,cap_term,np1,np2,counter, &
            max_ne, max_elem_indx,umb_node_counter,node_no,no_umb_art_nodes, umb_elem,&
            inlet_counter,outlet_counter,j,umb_elems_count,umb_elem_indx,ierror,ntemp,dwn_elems, &
            inlet1_np1,inlet2_np1,num_umb_nodes,umb_outlet_counter,inlet_with_gt_orders  
    integer, allocatable :: np_map(:)
    integer, allocatable :: ne_map(:)
    integer :: umb_art_outlets(2) = 0
    integer :: umb_art_outlet_nodes(2) = 0
    integer :: umb_ven_elems(3) = 0 !umbilical venous elements
    integer :: umb_ven_nodes(4) = 0 !umbilical venous nodes
    integer :: inlet_sord(2) = 0 !Strahler orders for arterial inlets
    integer :: inlet_hord(2) = 0 !Horsfield orders for arterial inlets    
    integer, allocatable :: umb_art_nodes(:)
    integer, allocatable :: umb_art_nodes_tmp(:)
    integer, allocatable :: umb_art_nodes_tmp2(:)
    integer :: umbilical_elems(20) = 0
    integer :: index, i
    character(LEN=132) :: ctemp1

    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'add_matching_mesh'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)
    
    if((umbilical_elem_option_in.NE.'same_as_arterial').AND.(umbilical_elem_option_in.NE.'single_umbilical_vein'))then
      print *, "Unknown option for creating a venous mesh. Using same_as_arterial."
      umbilical_elem_option = 'same_as_arterial'
    else
      umbilical_elem_option = umbilical_elem_option_in
    endif

    !if((umbilical_elem_option_in.EQ.'single_umbilical_vein').AND.(len(UMB_ELEMS_FILE).GT.0))then
    if((umbilical_elem_option_in.EQ.'single_umbilical_vein').AND.(size(umbilical_element_numbers).GT.0))then

       do index = 1, size(umbilical_element_numbers)
          umbilical_elems(index) = umbilical_element_numbers(index)
       enddo
       if(diagnostics_level.GT.1)then
          print *, "umbilical_elems", umbilical_elems
       endif

    endif

    umb_elems_count = count(umbilical_elems.NE.0)

    if((umbilical_elem_option_in.EQ.'single_umbilical_vein').AND.(umb_elems_count.EQ.0))then
      print *,"No umbilical elements found. Using option same_as_arterial."
      umbilical_elem_option = 'same_as_arterial'
    endif

    !Ultimately offset should be an input argument
    offset(1)=0.0_dp
    offset(2)=1e-6_dp
    offset(3)=0.0_dp

    if(umbilical_elem_option.EQ.'same_as_arterial')then

       ! the number of nodes after adding mesh will be:
       num_nodes_new = 2*num_nodes
       ! the number of elems after adding mesh will be:
       num_elems_new = 2*num_elems + num_units

    elseif(umbilical_elem_option.EQ.'single_umbilical_vein')then
       !populate the umbilical arterial nodes array and
       !find the inlet(s) and outlets in the umbilical elements

       !Count the number of nodes attached to the umbilical arteries
       allocate(umb_art_nodes_tmp(umb_elems_count*2))
       allocate(umb_art_nodes_tmp2(umb_elems_count*2))
       umb_node_counter = 0
       do noelem=1,umb_elems_count
         ne = umbilical_elems(noelem)
         do node_no=1,2
           umb_node_counter = umb_node_counter + 1
           np = elem_nodes(node_no,ne)
           umb_art_nodes_tmp(umb_node_counter) = np
         enddo
       enddo

       !The following 14 lines are modified from https://www.rosettacode.org/wiki/Remove_duplicate_elements
       umb_node_counter = 1
       umb_art_nodes_tmp2(1) = umb_art_nodes_tmp(1)
       outer: do i=2,size(umb_art_nodes_tmp)
         do j=1,umb_node_counter
          if (umb_art_nodes_tmp2(j) == umb_art_nodes_tmp(i)) then
           ! Found a match so start looking again
           cycle outer
        endif
        enddo
        ! No match found so add it to the output
        umb_node_counter = umb_node_counter + 1
        umb_art_nodes_tmp2(umb_node_counter) = umb_art_nodes_tmp(i)
       end do outer


       allocate(umb_art_nodes(umb_node_counter)) !ARC: we might be able to just define umb_art_nodes here from the above, but I leave as i havent checked if ordering matters
       deallocate(umb_art_nodes_tmp)
       deallocate(umb_art_nodes_tmp2)
       umb_art_nodes = 0

       inlet_counter = 0
       outlet_counter = 0
       umb_node_counter = 0
       do noelem=1,umb_elems_count
	  ne = umbilical_elems(noelem)
          !look for umbilical outlets = elements whose downstream elements are not umbilical elements themselves
          dwn_elems = elem_cnct(1,0,ne) !number of downstream elements
	  if(dwn_elems.GT.0)then
	     if(ALL(umbilical_elems.NE.elem_cnct(1,1,ne)))then
                outlet_counter = outlet_counter + 1
	        umb_art_outlets(outlet_counter) = ne
	        !store second nodes for mapping between arterial and venous elements
                umb_art_outlet_nodes(outlet_counter) = elem_nodes(2,ne)
	     endif
          endif

          do node_no=1,2
             np = elem_nodes(node_no,ne)
             if(ALL(umb_art_nodes.NE.np))then
                umb_node_counter = umb_node_counter + 1
                umb_art_nodes(umb_node_counter) = np
             endif
          enddo
       enddo

       if(count(umb_art_outlets.NE.0).EQ.0)then
          print *,"arterial umbilical outlets not found"
	  call exit(1)
       endif
       
       ! the number of nodes after adding mesh will be:
       num_nodes_new = 2*num_nodes - (umb_node_counter - 4) !umb_elems_count + 1 = number of nodes in umbilical section 
       !arteries = number of elems + 1; 4 umbilical vein nodes
       ! the number of elems after adding mesh will be:
       num_elems_new = 2*num_elems + num_units - (umb_elems_count - 3) !there are 3 umbilical venous elements
                                              
    endif

    allocate(np_map(num_nodes))
    np_map = 0
    allocate(ne_map(num_elems))
    ne_map = 0
    allocate(art_ven_elem_map(num_elems))
    !!! increase the size of node and element arrays to accommodate the additional elements
 
    call reallocate_node_elem_arrays(num_elems_new,num_nodes_new)
    ne0 = num_elems ! the starting local element number
    ne_global = num_elems ! assumes this is the highest element number (!!!)
    np0 = num_nodes ! the starting local node number
    np_global = num_nodes ! assumes this is the highest node number (!!!)

    umbilical_outlets = 0    
    if(umbilical_elem_option.EQ.'single_umbilical_vein')then

       node_counter = 1

       !a new venous node will be added at the centre of mass between copies of the last umbilical 
       !artery outlet nodes
       !if there are two inlets in the umbilical arterial elements, a venous node will be created 
       !at the centre of mass between umbilical arterial inlets + offset

       !create the umbilical vein nodes matching the umbilical arterial outlet nodes
       do indx=1,2
          np=np_global+node_counter
          umb_ven_nodes(indx) = np
          nonode = umb_art_outlet_nodes(indx)
          np_m=nodes(nonode)
          np_map(np_m)=np !maps new to old node numbering
          nodes(np0+node_counter)=np
          do nj=1,3
             node_xyz(nj,np)=node_xyz(nj,np_m)+offset(nj)
          enddo
          node_counter = node_counter + 1  
       enddo

       !create a node in between the venous outlet nodes
       !get the centre of mass between the two nodes
       np=np_global+node_counter
       umb_ven_nodes(3) = np
       nodes(np0+node_counter)=np
       do nj=1,3
          node_xyz(nj,np)=(node_xyz(nj,umb_ven_nodes(1))+node_xyz(nj,umb_ven_nodes(2)))/2
       enddo   
       node_counter = node_counter + 1  

       !if two inlets then create a new node in between them
       if(count(umbilical_inlets.NE.0).GT.1)then
          np=np_global+node_counter
          umb_ven_nodes(4) = np
          nodes(np0+node_counter)=np
          inlet1_np1 = elem_nodes(1,umbilical_inlets(1)) 
          inlet2_np1 = elem_nodes(1,umbilical_inlets(2))    
          do nj=1,3
             node_xyz(nj,np)=(node_xyz(nj,inlet1_np1)+node_xyz(nj,inlet2_np1))/2 + offset(nj)
          enddo
          node_counter = node_counter + 1  

       else
          !if one inlet copy the first node
          np=np_global+node_counter
          umb_ven_nodes(4) = np
          nodes(np0+node_counter)=np  
          np1 = elem_nodes(1,umbilical_inlets(1)) 
          do nj=1,3
             node_xyz(nj,np)=node_xyz(nj,np1) + offset(nj)
          enddo
          node_counter = node_counter + 1  

       endif

       !copy non umbilical arterial nodes
       do np_m=1,num_nodes
          if(ALL(umb_art_nodes.NE.np_m))then
             np=np_global+node_counter
             !np_m=nodes(nonode)
             np_map(np_m)=np !maps new to old node numbering
             nodes(np0+node_counter)=np
             do nj=1,3
               node_xyz(nj,np)=node_xyz(nj,np_m)+offset(nj)
             enddo
             node_counter = node_counter + 1
          endif
       enddo
    else !umbilical_elem_option = same_as_arterial

       do np_m=1,num_nodes
          np=np_global+np_m
          !np_m=nodes(nonode)
          np_map(np_m)=np !maps new to old node numbering
          nodes(np0+np_m)=np
          do nj=1,3
            node_xyz(nj,np)=node_xyz(nj,np_m)+offset(nj)
          enddo
        !Doesnt map versions, would be added here
       enddo

    endif !umbilical_elem_option

    elems_at_node(np,0)=0 !initialise

    elem_counter = 1

    if(umbilical_elem_option.EQ.'single_umbilical_vein')then

       !create an element between nodes 3 and 4
       !umbilical venous outlet element
       ne=ne_global+elem_counter
       umb_ven_elems(1)=ne
       umbilical_outlets(1) = ne !store umbilical venous outlet
       elem_field(ne_group,ne)=2.0_dp!VEIN
       elems(ne0+elem_counter)=ne
       elem_nodes(1,ne)=umb_ven_nodes(3)
       elem_nodes(2,ne)=umb_ven_nodes(4)
       elem_counter = elem_counter + 1      

       !create an element between nodes 1 and 3
       ne=ne_global+elem_counter
       umb_ven_elems(2)=ne
       elem_field(ne_group,ne)=2.0_dp!VEIN
       elems(ne0+elem_counter)=ne
       elem_nodes(1,ne)=umb_ven_nodes(1)
       elem_nodes(2,ne)=umb_ven_nodes(3)
       elem_counter = elem_counter + 1

       !create an element between nodes 2 and 3
       ne=ne_global+elem_counter
       umb_ven_elems(3)=ne
       elem_field(ne_group,ne)=2.0_dp!VEIN
       elems(ne0+elem_counter)=ne
       elem_nodes(1,ne)=umb_ven_nodes(2)
       elem_nodes(2,ne)=umb_ven_nodes(3)
       elem_counter = elem_counter + 1
     

       !populate elements at a node for the new umbilical venous elements
       do indx=1,3
          noelem=umb_ven_elems(indx)
          elems_at_node(elem_nodes(1,noelem),0)=elems_at_node(elem_nodes(1,noelem),0)+1
          elems_at_node(elem_nodes(1,noelem),elems_at_node(elem_nodes(1,noelem),0))=noelem
          elems_at_node(elem_nodes(2,noelem),0)=elems_at_node(elem_nodes(2,noelem),0)+1
          elems_at_node(elem_nodes(2,noelem),elems_at_node(elem_nodes(2,noelem),0))=noelem
        enddo
    endif

    umb_outlet_counter = 0
    do ne_m=1,num_elems
        !ne_m=elems(noelem)
        elem_field(ne_group,ne_m)=0.0_dp!ARTERY
        if((umbilical_elem_option.EQ.'same_as_arterial').OR. &
                    ((umbilical_elem_option.EQ.'single_umbilical_vein').AND.(ALL(umbilical_elems.NE.ne_m))))then
           ne=ne_global+elem_counter
           !if copying arterial umbilical inlets, store venous umbilical outlets
           if(.not.ALL(umbilical_inlets.NE.ne_m))then
              umb_outlet_counter = umb_outlet_counter + 1
              umbilical_outlets(umb_outlet_counter) = ne
           endif
           elem_field(ne_group,ne)=2.0_dp!VEIN
           elems(ne0+elem_counter)=ne
           ne_map(ne_m)=ne !maps new to old element numbering
           elem_nodes(1,ne)=np_map(elem_nodes(2,ne_m))
           elem_nodes(2,ne)=np_map(elem_nodes(1,ne_m))   
           !if worrying about regions and versions do it here
           elems_at_node(elem_nodes(1,ne),0)=elems_at_node(elem_nodes(1,ne),0)+1
           elems_at_node(elem_nodes(1,ne),elems_at_node(elem_nodes(1,ne),0))=ne
           elems_at_node(elem_nodes(2,ne),0)=elems_at_node(elem_nodes(2,ne),0)+1
           elems_at_node(elem_nodes(2,ne),elems_at_node(elem_nodes(2,ne),0))=ne
           
           nindex=no_gen
           elem_ordrs(nindex,ne)=elem_ordrs(nindex,ne_m)
           nindex=no_sord
           elem_ordrs(nindex,ne)=elem_ordrs(nindex,ne_m)
           nindex=no_hord
           elem_ordrs(nindex,ne)=elem_ordrs(nindex,ne_m)

           elem_counter = elem_counter + 1
        endif
    enddo 

    !populate mapping between arterial and venous elements (for radius from file option)
    art_ven_elem_map = ne_map

    if(diagnostics_level.GE.1)then
       do umb_outlet_counter=1,2
          if(umbilical_outlets(umb_outlet_counter).GT.0)then
             print *,"venous umbilical outlet element number",umbilical_outlets(umb_outlet_counter)
          endif
       enddo
    endif
  
    max_ne = ne
    max_elem_indx = ne0 + elem_counter - 1

    !copy element connectivity
    do ne_m=1,num_elems
           !ne_m=elems(noelem)
           ne=ne_map(ne_m)
           if(ne.GT.0)then             
              elem_cnct(-1,0,ne)=elem_cnct(1,0,ne_m)
              elem_cnct(1,0,ne)=elem_cnct(-1,0,ne_m)
              do n=1,elem_cnct(-1,0,ne)     
                elem_cnct(-1,n,ne)=ne_map(elem_cnct(1,n,ne_m))
              enddo
              do n=1,elem_cnct(1,0,ne) 
                elem_cnct(1,n,ne)=ne_map(elem_cnct(-1,n,ne_m))
              enddo    
           endif   
    enddo

    if(umbilical_elem_option.EQ.'single_umbilical_vein')then

        !populate element connectivity for the new umbilical venous elements
        elem_cnct(-1,0,umb_ven_elems(1)) = 2
        elem_cnct(-1,1,umb_ven_elems(1)) = umb_ven_elems(2)
        elem_cnct(-1,2,umb_ven_elems(1)) = umb_ven_elems(3)
        elem_cnct(1,0,umb_ven_elems(1)) = 0 

        elem_cnct(-1,0,umb_ven_elems(2)) = 0
        elem_cnct(1,0,umb_ven_elems(2)) = 1
        elem_cnct(1,1,umb_ven_elems(2)) = umb_ven_elems(1)

        elem_cnct(-1,0,umb_ven_elems(3)) = 0
        elem_cnct(1,0,umb_ven_elems(3)) = 1
        elem_cnct(1,1,umb_ven_elems(3)) = umb_ven_elems(1)


        !fix element connectivity for elements upstream of new umbilical venous elements  
        do indx=2,3
           noelem = umb_ven_elems(indx)
           nonode = elem_nodes(1,noelem) !get the first node
	   !get elements connected to the first node
           do indx2=1,elems_at_node(nonode,0)
	      elem_at_node = elems_at_node(nonode,indx2)             
              if(elem_at_node.NE.noelem)then
          
                 elem_cnct(1,0,elem_at_node)=1
	         elem_cnct(1,1,elem_at_node)=noelem

                 elem_cnct(-1,0,noelem)=elem_cnct(-1,0,noelem)+1
                 elem_cnct(-1,elem_cnct(-1,0,noelem),noelem)=elem_at_node
              endif	
	   enddo
        enddo

        !element orders for the first new element - the same as the arterial inlet
        !if two inlets, copy greater Strahler and Horsfield orders; generation for 
        !each inlet is 1
        inlet_with_gt_orders = 1
        if(count(umbilical_inlets.NE.0).GT.1)then
           do indx=1,2
              inlet_sord(indx) = elem_ordrs(no_sord,umbilical_inlets(indx))
              inlet_hord(indx) = elem_ordrs(no_hord,umbilical_inlets(indx))
           enddo
           if((inlet_sord(1).LT.inlet_sord(2)).OR.(inlet_hord(1).LT.inlet_hord(2)))then
              inlet_with_gt_orders = 2
           endif
        endif
        nindex=no_gen
        elem_ordrs(nindex,umb_ven_elems(1))=elem_ordrs(nindex,umbilical_inlets(inlet_with_gt_orders))
        nindex=no_sord
        elem_ordrs(nindex,umb_ven_elems(1))=elem_ordrs(nindex,umbilical_inlets(inlet_with_gt_orders))
        nindex=no_hord
        elem_ordrs(nindex,umb_ven_elems(1))=elem_ordrs(nindex,umbilical_inlets(inlet_with_gt_orders))

        !populate mapping between arterial and venous elements (for radius from file option)
        art_ven_elem_map(umbilical_inlets(inlet_with_gt_orders)) = umb_ven_elems(1)

        !element orders for the second and third new element - the same as the arterial 
        !umbilical outlets
        do indx=2,3
           nindex=no_gen
           elem_ordrs(nindex,umb_ven_elems(indx))=elem_ordrs(nindex,umb_art_outlets(indx-1))
           nindex=no_sord
           elem_ordrs(nindex,umb_ven_elems(indx))=elem_ordrs(nindex,umb_art_outlets(indx-1))
           nindex=no_hord
           elem_ordrs(nindex,umb_ven_elems(indx))=elem_ordrs(nindex,umb_art_outlets(indx-1))

           !populate mapping between arterial and venous elements (for radius from file option)
           art_ven_elem_map(umb_art_outlets(indx-1)) = umb_ven_elems(indx)
        enddo

     endif !umbilical_elem_option

     cap_conns=0
     cap_term=0
     noelem0 = max_elem_indx
     ne1 = max_ne
     do nu=1,num_units
         ne=units(nu)
         cap_term=cap_term+1
         np1=elem_nodes(2,ne)
         np2=np_map(np1)
         noelem0=noelem0+1
         ne1=ne1+1
         elems(noelem0)=ne1
         elem_nodes(1,ne1)=np1
         elem_nodes(2,ne1)=np2
         elems_at_node(np1,0)=elems_at_node(np1,0)+1
         elems_at_node(np1,elems_at_node(np1,0))=ne1
         elems_at_node(np2,0)=elems_at_node(np2,0)+1
         elems_at_node(np2,elems_at_node(np2,0))=ne1
         elem_cnct(1,0,ne)=elem_cnct(1,0,ne)+1 
         elem_cnct(1,elem_cnct(1,0,ne),ne)=ne1       
         ne_m=ne
         elem_cnct(-1,0,ne_map(ne_m))=elem_cnct(-1,0,ne_map(ne_m))+1
         elem_cnct(-1,elem_cnct(-1,0,ne_map(ne_m)),ne_map(ne_m))=ne1
         elem_cnct(-1,0,ne1)=1
         elem_cnct(1,0,ne1)=1
         elem_cnct(-1,1,ne1)=ne
         elem_cnct(1,1,ne1)=ne_map(ne_m)
         nindex=no_gen
         elem_ordrs(nindex,ne1)=elem_ordrs(nindex,ne_m)
         nindex=no_sord
         elem_ordrs(nindex,ne1)=elem_ordrs(nindex,ne_m)
         nindex=no_hord
         elem_ordrs(nindex,ne1)=elem_ordrs(nindex,ne_m)
         elem_field(ne_group,ne1)=1.0_dp!connection between meshes
    enddo
    if(diagnostics_level.GT.1)then
      print *, 'Number of connections', cap_term
    endif
 
    num_nodes=num_nodes_new
    num_arterial_elems = num_elems
    num_elems=num_elems_new
  
    !calculate the element lengths and directions for the venous elements
    do ne=num_arterial_elems+1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       elem_field(ne_length,ne) = DSQRT((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)   
       do j=1,3
          elem_direction(j,ne) = (node_xyz(j,np2) - &
               node_xyz(j,np1))/elem_field(ne_length,ne)           
       enddo

    enddo

  
    if(diagnostics_level.GT.1)then 
        print *, "num_nodes=",num_nodes
        print *, "num_arterial_elems=",num_arterial_elems
        print *, "num_elems=",num_elems
   	 	!print out new node geometry, number of connected elements for each element and element connectivity array
   		print *,"element nodes:"
   		DO ne=1,num_elems
   	   		DO nj=1,2 !each element has 2 nodes in 1D tree
   	     		print *, "elem_nodes(nj,ne)",nj,ne,"=",elem_nodes(nj,ne)
   	   		ENDDO
   		ENDDO
   
    		DO n=1,num_nodes
       		DO nj=1,3
         		print *,"node_xyz(nj,n)",nj,n,"=",node_xyz(nj,n)
       		ENDDO
    		ENDDO
   
    		DO n=1,num_nodes
       		print *," "
       		print *,"node",n
       		print *, "total number of elements connected",elems_at_node(n,0)
       		DO ne=1,elems_at_node(n,0)
          		print *, "element",elems_at_node(n,ne)
       		ENDDO
    		ENDDO
 
		! total count of upstream elements connected to element ne elem_cnct(-1,0,ne)
		! upstream elements elem_cnct(-1,counter,ne)
		! total count of downstream elements connected to element ne elem_cnct(1,0,ne)
		! downstream elements elem_cnct(1,counter,ne)
   		DO ne=1,num_elems
   	    		print *,""
   	    		print *,"element",ne 
       		IF(elem_cnct(-1,0,ne).gt.0)THEN
       	    		print *, "total number of upstream elements:",elem_cnct(-1,0,ne)
       			DO counter=1,elem_cnct(-1,0,ne)
          			print *, "upstream element",elem_cnct(-1,counter,ne)
       	    		ENDDO
       		ENDIF
       		IF(elem_cnct(1,0,ne).gt.0)THEN
       	    		print *, "total number of downstream elements:",elem_cnct(1,0,ne)
       			DO counter=1,elem_cnct(1,0,ne)
          			print *, "downstream element",elem_cnct(1,counter,ne)
       	    		ENDDO
       		ENDIF    
    		ENDDO   
    		
    		DO ne=1,num_elems
    			do counter=1,3		
    		  		print *, "elem_ordrs(",counter,",",ne,")=",elem_ordrs(counter,ne)
    			ENDDO
    		enddo
    		 
    endif !diagnostics_level
       
    deallocate(np_map)
    call enter_exit(sub_name,2)

  end subroutine add_matching_mesh

!
!###################################################################################
!
  subroutine append_units()
  !*Description:* Defines terminal units at the end of a tree structure
    use arrays,only: dp, elem_cnct,elem_symmetry,elem_units_below,&
         num_elems,num_units,units,unit_field
    use indices,only: num_nu
    use diagnostics, only: enter_exit,get_diagnostics_level
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_APPEND_UNITS" :: APPEND_UNITS

    integer :: ne,ne0,nu
    character(len=60) :: sub_name
    integer:: diagnostics_level

    sub_name = 'append_units'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    num_units = 0
    DO ne=1,num_elems
       IF(elem_cnct(1,0,ne).eq.0)THEN! terminal element
          num_units=num_units+1
       ENDIF
    ENDDO
    
    if(diagnostics_level.GT.1)then
		print *,"num_units=",num_units
	endif

    if(allocated(units))then !increasing the array size; just overwrite
       deallocate(units)
       deallocate(unit_field)
    endif
    allocate(units(num_units))
    allocate(unit_field(num_nu,num_units))

    unit_field=0.0_dp
    units=0
    elem_units_below(1:num_elems) = 0 !initialise the number of terminal units below a branch

    nu=0
    DO ne=1,num_elems
       IF(elem_cnct(1,0,ne).eq.0)THEN
          nu=nu+1
          units(nu)=ne     !Set up units array containing terminals
          elem_units_below(ne)=1
       ENDIF
    ENDDO

    ! count the effective number of elements below each branch
    do ne=num_elems,1,-1
       ne0=elem_cnct(-1,1,ne)
       if(ne0.gt.0)then !not an inlet
       elem_units_below(ne0) = elem_units_below(ne0) &
            + elem_units_below(ne)*elem_symmetry(ne)
       endif
    enddo !ne
	
	if(diagnostics_level.GT.1)then
		do ne=1,num_elems
			print *,"elem_units_below(",ne,")",elem_units_below(ne)
		enddo
	endif

    call enter_exit(sub_name,2)

  end subroutine append_units

!
!###################################################################################
!
  subroutine calc_capillary_unit_length(num_convolutes,num_generations)
  !*Description:* Calculates the effective length of a capillary unit based on its total resistance
  ! and assumed radius, given the number of terminal convolute connections and the number of
  ! generations of symmetric intermediate villous branches 
    use arrays,only: dp,num_units,units,elem_field,elem_direction, &
                     node_xyz,elem_nodes,elem_cnct,num_conv,num_conv_gen, &
                     cap_resistance,terminal_resistance,terminal_length, &
                     cap_radius,is_capillary_unit,total_cap_volume,total_cap_surface_area, &
                     num_elems
    use diagnostics, only: enter_exit,get_diagnostics_level
    use other_consts, only: PI
    use indices, only: ne_length,ne_radius,ne_vol
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CALC_CAPILLARY_UNIT_LENGTH" :: CALC_CAPILLARY_UNIT_LENGTH

    integer, intent(inout) :: num_convolutes,num_generations

    real(dp) :: int_length,int_radius,seg_length,viscosity, &
                seg_resistance,cap_unit_radius, cap_length, &
                capillary_unit_vol, capillary_unit_area
    real(dp) :: int_radius_in,int_radius_out,avg_art_rad,avg_ven_rad
    real(dp),allocatable :: resistance(:)
    integer :: ne,nu,i,j,np1,np2,nc,nv
    integer :: AllocateStatus
    character(len=60) :: sub_name
    integer:: diagnostics_level

    sub_name = 'calc_terminal_unit_length'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    !check number of capillary convolutes and number of intermediate villous tree generations
    if (num_convolutes.LE.0)then
      num_convolutes = 6
    endif
    if (num_generations.LE.0)then
      num_generations = 3
    endif

    num_conv = num_convolutes
    num_conv_gen = num_generations

    if (diagnostics_level.GE.2)then
      print *, "num_convolutes=",num_convolutes
      print *, "num_generations=",num_generations
    endif

    allocate (resistance(num_convolutes+1), STAT = AllocateStatus)
    if (AllocateStatus /= 0)then
       STOP "*** Not enough memory for resistance array ***"
    endif

    allocate(is_capillary_unit(num_elems), STAT = AllocateStatus)
    if (AllocateStatus /= 0)then
       STOP "*** Not enough memory for is_capillary_unit array ***"
    endif

    avg_art_rad = 0.0_dp
    avg_ven_rad = 0.0_dp
    do nu=1,num_units
       ne =units(nu) !Get a terminal unit
       nc = elem_cnct(1,1,ne) !capillary unit is downstream of a terminal unit
       nv =  elem_cnct(1,1,nc) !vein is downstream of the capillary
       avg_art_rad = avg_art_rad + elem_field(ne_radius,ne)
       avg_ven_rad = avg_ven_rad + elem_field(ne_radius,nv)
    enddo
    avg_art_rad = avg_art_rad/(num_units*1.0_dp)
    avg_ven_rad = avg_ven_rad/(num_units*1.0_dp)
    int_radius_in = (avg_art_rad+avg_ven_rad)/2.0_dp ! mm radius of inlet intermediate villous (average of artery and vein)
    int_radius_out=(0.03_dp + 0.03_dp/2.0_dp)/2.0_dp ! mm radius of mature intermediate villous (average of artery and vein)
    int_length=1.5_dp !mm Length of each intermediate villous
    cap_length=3_dp/num_convolutes !mm length of capillary convolutes
    cap_radius=0.0144_dp/2.0_dp !radius of capillary convolutes
    seg_length=int_length/num_convolutes !lengh of each intermediate villous segment
    viscosity=0.33600e-02_dp !Pa.s !viscosity: fluid viscosity
    cap_unit_radius = 0.03_dp
    cap_resistance=(8.d0*viscosity*cap_length)/(PI*cap_radius**4) !resistance of each capillary convolute segment
    terminal_resistance = 0

    !calculate total capillary unit volume and surface area so that this can be used by subroutine calculate_stats
    !in module pressure_resistance_flow
    capillary_unit_vol = PI * cap_radius**2 * cap_length * num_convolutes * num_generations
    capillary_unit_area = 2 * PI * cap_radius * cap_length * num_convolutes * num_generations

    do j=1,num_generations

      int_radius = int_radius_in - (int_radius_in-int_radius_out)/num_generations*j

      !capillary unit volume and surface area calculation - adding intermediate villous volume 
      !and area for the generation
      capillary_unit_vol = capillary_unit_vol + PI * int_radius**2 * int_length
      capillary_unit_area = capillary_unit_area + 2 * PI * int_radius * int_length

      seg_resistance=(8.d0*viscosity*seg_length)/(PI*int_radius**4) !resistance of each intermediate villous segment
      !calculate total resistance of terminal capillary conduits
      i=1
      resistance(i)= cap_resistance + 2.d0*seg_resistance
      do i=2,num_convolutes
        resistance(i)=2.d0*seg_resistance + 1/(1/cap_resistance + 1/resistance(i-1))
      enddo
      cap_resistance = resistance(num_convolutes) !Pa . s per mm^3 total resistance of terminal capillary conduits
	
      !We have symmetric generations of intermediate villous trees so we can calculate the total resistance
      !of the system by summing the resistance of each generation

      terminal_resistance = terminal_resistance + cap_resistance/2**j
    enddo

    terminal_length = terminal_resistance*(PI*cap_unit_radius**4)/(8.d0*viscosity)

    total_cap_volume = capillary_unit_vol * num_units/1000 !in cm**3
    total_cap_surface_area = capillary_unit_area * num_units/100 !in cm**2

    if(diagnostics_level.GE.2)then
      print *, "Resistance of capillary conduits=",cap_resistance
      print *, "Resistance of all generations of capillaries per terminal unit=",terminal_resistance
      print *, "Effective length of each capillary unit (mm)",terminal_length
      print *, "Total capillary length for the vasculature (cm)",(terminal_length*num_units)/10
      print *, "Total capillary volume (cm**3) = ",total_cap_volume
      print *, "Total capillary surface area (cm**2) = ", total_cap_surface_area
    endif

    is_capillary_unit = 0
    !set the effective length of each capillary unit based the total resistance of capillary convolutes   
    do nu=1,num_units
      ne =units(nu) !Get a terminal unit    
      nc = elem_cnct(1,1,ne) !capillary unit is downstream of a terminal unit
      is_capillary_unit(nc) = 1
      !update element radius
      elem_field(ne_radius,nc) = cap_unit_radius
      !update element length   
      elem_field(ne_length,nc) = terminal_length
      !update element direction
      np1=elem_nodes(1,nc)
      np2=elem_nodes(2,nc)
      do j=1,3
        elem_direction(j,nc) = (node_xyz(j,np2) - &
               node_xyz(j,np1))/elem_field(ne_length,nc)           
      enddo

      !update element volume
      elem_field(ne_vol,nc) = PI * elem_field(ne_radius,nc)**2 * &
            elem_field(ne_length,nc)

    enddo

    call enter_exit(sub_name,2)

  end subroutine calc_capillary_unit_length
!
!###################################################################################
!
  subroutine define_1d_elements(ELEMFILE,anastomosis_elem_in)
  !*Description:* Reads in an element ipelem file to define a geometry
    use arrays,only: dp, elem_direction,elem_field,elems,elem_cnct,elem_nodes,&
         elem_ordrs,elem_symmetry,elems_at_node,elem_units_below,&
         node_xyz,num_elems,num_nodes,elem_cnct_no_anast,anastomosis_elem, &
         umbilical_inlets
    use indices
    use diagnostics, only: enter_exit,get_diagnostics_level
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_1D_ELEMENTS" :: DEFINE_1D_ELEMENTS

    character(len=MAX_FILENAME_LEN), intent(in) :: ELEMFILE
    integer, optional :: anastomosis_elem_in
    !     Local Variables
    integer :: ibeg,iend,ierror,i_ss_end,j,ne,ne_global,&
         nn,np,np1,np2,np_global,inlet_counter
    character(LEN=132) :: ctemp1
    character(LEN=40) :: sub_string
    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'define_1d_elements'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    anastomosis_elem = anastomosis_elem_in

    open(10, file=ELEMFILE, status='old')

    read_number_of_elements : do
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "elements")> 0) then
          call get_final_integer(ctemp1,num_elems)
          if(diagnostics_level.GT.1)then
          	print *, "num_elems", num_elems
          endif
          exit read_number_of_elements
       endif
    end do read_number_of_elements

!!! allocate memory for element arrays
    if(allocated(elems)) deallocate(elems)
    allocate(elems(num_elems))
    if(allocated(elem_cnct)) deallocate(elem_cnct)
    allocate(elem_cnct(-1:1,0:2,0:num_elems))
    if(allocated(elem_cnct_no_anast)) deallocate(elem_cnct_no_anast)
    allocate(elem_cnct_no_anast(-1:1,0:2,0:num_elems))
    if(allocated(elem_nodes)) deallocate(elem_nodes)
    allocate(elem_nodes(2,num_elems))
    if(allocated(elem_ordrs)) deallocate(elem_ordrs)
    allocate(elem_ordrs(num_ord,num_elems))
    if(allocated(elem_symmetry)) deallocate(elem_symmetry)
    allocate(elem_symmetry(num_elems))
    if(allocated(elem_units_below)) deallocate(elem_units_below)
    allocate(elem_units_below(num_elems))
    if(allocated(elems_at_node)) deallocate(elems_at_node)
    allocate(elems_at_node(num_nodes,0:3))
    if(allocated(elem_field)) deallocate(elem_field)
    allocate(elem_field(num_ne,num_elems))
    if(allocated(elem_direction)) deallocate(elem_direction)
    allocate(elem_direction(3,num_elems))

!!! initialise element arrays
    elems=0
    elem_nodes=0
    elem_symmetry = 1
    elem_field = 0.0_dp

    ne=0
    !each element has 2 nodes in 1D tree - this is defined in array elem_nodes: elem_nodes (1,element_number ne) = node np"	
    read_an_element : do
       !.......read element number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Element")> 0) then
          call get_final_integer(ctemp1,ne_global) !get element number
          ne=ne+1
          elems(ne)=ne_global
             read_element_nodes : do
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "global")> 0) then !found the correct line
                iend=len(ctemp1)
                ibeg=index(ctemp1,":")+1 !get location of first integer in string
                sub_string = adjustl(ctemp1(ibeg:iend)) ! get the characters beyond : remove leading blanks
                i_ss_end=len(sub_string) !get the end location of the sub-string
                ibeg=1
                do nn=1,2
                   iend=index(sub_string," ") !get location of first blank in sub-string
                   read (sub_string(ibeg:iend-1), '(i7)' ) np_global
                   call get_local_node(np_global,np) ! get local node np for global node
                   elem_nodes(nn,ne)=np ! the local node number, not global             
                   if(diagnostics_level.GT.1)then
                   		print *,"elem_nodes(nn,ne)", nn, ne, "= np", np
                   endif
                   sub_string = adjustl(sub_string(iend:i_ss_end)) ! get chars beyond blank, remove leading blanks
                end do
                exit read_element_nodes
             endif !index
          end do read_element_nodes
          if(ne.ge.num_elems) exit read_an_element
       endif

    end do read_an_element

    close(10)

    ! calculate the element lengths and directions
    do ne=1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       elem_field(ne_length,ne) = DSQRT((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)   
       do j=1,3
          elem_direction(j,ne) = (node_xyz(j,np2) - &
               node_xyz(j,np1))/elem_field(ne_length,ne)           
       enddo !j

    enddo

    call element_connectivity_1d


    !populate umbilical_inlets array
    inlet_counter = 0
    umbilical_inlets = 0
    do ne=1,num_elems
       if(elem_cnct(-1,0,ne).EQ.0)then
          inlet_counter = inlet_counter + 1
          umbilical_inlets(inlet_counter) = ne
          if(diagnostics_level.GE.1)then
             print *,"umbilical inlet element number =",ne
          endif
       endif
    enddo
    if(count(umbilical_inlets.NE.0).EQ.0)then
       print *,"inlet not found"
       call exit(1)
    endif

    call evaluate_ordering

    call enter_exit(sub_name,2)

  END subroutine define_1d_elements

!
!###################################################################################
!
  subroutine define_node_geometry(NODEFILE)
  !*Description:* Reads in an ipnode file to define a tree geometry
    use arrays,only: dp,nodes,node_field,node_xyz,num_nodes
    use diagnostics, only: enter_exit,get_diagnostics_level
    use indices
    use other_consts, only: MAX_FILENAME_LEN
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_NODE_GEOMETRY" :: DEFINE_NODE_GEOMETRY

    character(len=MAX_FILENAME_LEN), intent(in) :: NODEFILE !Input nodefile
    !     Local Variables
    integer :: i,ierror,np,np_global,&
         num_versions,nv,NJT
    character(LEN=132) :: ctemp1
    LOGICAL :: versions
    real(dp) :: point
    character(len=60) :: sub_name
    integer :: diagnostics_level
  
    sub_name = 'define_node_geometry'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)
   

    versions = .TRUE.
    NJT = 0
    open(10, file=NODEFILE, status='old')

    !.....read in the total number of nodes. read each line until one is found
    !.....that has the correct keyword (nodes). then return the integer that is
    !.....at the end of the line
    read_number_of_nodes : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "nodes")> 0) then !keyword "nodes" is found in ctemp1
          call get_final_integer(ctemp1,num_nodes) !return the final integer
          if(diagnostics_level.GT.1)then
       	  	print *, "num_nodes",num_nodes
       	  endif
          exit read_number_of_nodes !exit the named do loop
       endif

    end do read_number_of_nodes

    if(allocated(nodes)) deallocate (nodes)
    allocate (nodes(num_nodes))
    if(allocated(node_xyz)) deallocate (node_xyz)
    allocate (node_xyz(3,num_nodes))
    if(allocated(node_field)) deallocate (node_field)
    allocate (node_field(num_nj,num_nodes))
    nodes = 0 !initialise node index values
    node_xyz = 0.0_dp !initialise
    node_field = 0.0_dp !initialise

    !.....read in the number of coordinates
    read_number_of_coords : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "coordinates")> 0) then !keyword "coordinates" is found
          call get_final_integer(ctemp1,NJT) !return the final integer
          exit read_number_of_coords !exit the named do loop
       endif
    end do read_number_of_coords

    !.....check whether versions are prompted (>1)
    read_versions : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "different")> 0) then !keyword "different" is found
          if(index(ctemp1, " N")> 0) then !keyword " N" is found
             versions=.false.
          endif
          exit read_versions !exit the named do loop
       endif
    end do read_versions

!!! WARNING :: following should be in general code
    ! note that only the first version of coordinate is currently read in

    !.....read the coordinate, derivative, and version information for each node.
    np=0
    read_a_node : do !define a do loop name
       !.......read node number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Node")> 0) then
          call get_final_integer(ctemp1,np_global) !get node number
          np=np+1
          nodes(np)=np_global
          !.......read coordinates and derivatives
          do i=1,NJT ! for the NJT coordinates
             !...........coordinate
             num_versions=1
             if(versions)then
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                call get_final_integer(ctemp1,num_versions)
             endif
             if(num_versions > 1)then
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !temporary line
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                call get_final_real(ctemp1,point)
                do nv=2,num_versions
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !temporary line
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                enddo
             else
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                call get_final_real(ctemp1,point)
             endif
             node_xyz(i,np)=point
             if(diagnostics_level.GT.1)then
             	print *, "node_xyz(i,np)",i," ",np,"=", point
             endif
          end do !i

       endif !index
                  
       if(np.ge.num_nodes) exit read_a_node
    end do read_a_node

    close(10)
    
    call enter_exit(sub_name,2)

  END subroutine define_node_geometry
!
!###################################################################################
!
  subroutine define_rad_from_file(FIELDFILE,order_system, s_ratio)
  !*Description:* Reads in a radius field associated with the full feto-placental tree
  ! and assigns radius information to each element. If only radii for part of the tree
  ! are available (starting with the inlet), the remaining arterial vessel radii are 
  ! calcuated based on the order system and ratio in the input arguments. If venous
  ! vessels were added using add_matching_mesh and their radius was not available in 
  ! the radius field file, define_rad_from_geometry must be used for the venous vessels.

    use arrays,only: dp,elem_field,elem_cnct,elem_nodes,&
         elems_at_node,num_elems,num_nodes, num_arterial_elems, elem_ordrs, &
         umbilical_inlets
    use indices,only: ne_length,ne_radius,ne_radius_in,ne_radius_out,no_sord, &
                    no_hord,ne_vol 
    use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
    use diagnostics, only: enter_exit,get_diagnostics_level
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_RAD_FROM_FILE" :: DEFINE_RAD_FROM_FILE

    character(len=MAX_FILENAME_LEN), intent(in) :: FIELDFILE
    character(len=MAX_STRING_LEN), optional ::  order_system
    real(dp), optional :: s_ratio

    !     Local Variables
    real(dp), allocatable :: node_radius(:)
    integer :: ierror,ne,np,np1,np2,np_global,surround,radii_num_nodes, &
               minelem_no_radius,start_elem,n_max_ord,nindex,remaining_elems, &
               indx
    character(LEN=132) :: ctemp1
    LOGICAL :: versions
    real(dp) :: radius,start_rad
    character(len=MAX_STRING_LEN) ::  venous_option
    integer :: inlet_ord(2) = 0 !Strahler or Horsfield orders for arterial inlets
    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'define_rad_from_file'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)
    
    versions = .TRUE.

    open(10, file=FIELDFILE, status='old')

    !.....read the number of nodes for which we have radii
    !.....read each line until one is found
    !.....that has the correct keyword (nodes). then return the integer that is
    !.....at the end of the line
    read_number_of_nodes : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "nodes")> 0) then !keyword "nodes" is found in ctemp1
          call get_final_integer(ctemp1,radii_num_nodes) !return the final integer
          if(diagnostics_level.GT.1)then
       	  	print *, "radii_num_nodes",radii_num_nodes
       	  endif
          exit read_number_of_nodes !exit the named do loop
       endif

    end do read_number_of_nodes


   allocate(node_radius(radii_num_nodes))

    !.....check whether versions are prompted (>1) - versions are not used later in this subroutine
    read_versions : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "different")> 0) then !keyword "different" is found
          if(index(ctemp1, " N")> 0) then !keyword " N" is found
             versions=.false.
          endif
          exit read_versions !exit the named do loop
       endif
    end do read_versions

    np = 0
    !.....read the coordinate, derivative, and version information for each node.
    read_a_node : do !define a do loop name
       !.......read node number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Node")> 0) then
          call get_final_integer(ctemp1,np_global) !get global node number
          ! find the corresponding local node number
          call get_local_node(np_global,np) ! get local node np for global node
          read(unit=10, fmt="(a)", iostat=ierror) ctemp1
          read(unit=10, fmt="(a)", iostat=ierror) ctemp1
    	  if(index(ctemp1, "value")> 0) then
                call get_final_real(ctemp1,radius)
                node_radius(np)=radius   
          endif
       endif !index
       if(np.ge.radii_num_nodes) exit read_a_node
    end do read_a_node

    minelem_no_radius = num_elems
    ! assign a radius to each element whose 2nd node is in the radius file
    do ne=1,num_elems
       !for each element set the element radius to the radius of the second node
       np2=elem_nodes(2,ne)
       if(np2.LE.radii_num_nodes)then
          elem_field(ne_radius,ne) = node_radius(np2)
          elem_field(ne_radius_in,ne) = node_radius(np2)
          elem_field(ne_radius_out,ne) = node_radius(np2)
          if(diagnostics_level.GT.1)then
     	     print *,"radius for element",ne,"=",elem_field(ne_radius,ne) 
     	     print *,"radius in for element",ne,"=",elem_field(ne_radius_in,ne)
     	     print *,"radius out for element",ne,"=",elem_field(ne_radius_out,ne)      
          endif
       elseif(minelem_no_radius.NE.0)then
          if(ne.LT.minelem_no_radius)then
             minelem_no_radius = ne
          endif
       endif
    enddo

    if(minelem_no_radius.LT.num_elems)then

       ! calculate radii for the remaining arterial elements whose 2nd node radius 
       ! was not in the file
       if(num_arterial_elems.GT.0)then
          remaining_elems = num_arterial_elems
       else
          remaining_elems = num_elems
       endif

       !Strahler and Horsfield ordering system
       if(ORDER_SYSTEM(1:5).EQ.'strah')THEN
          nindex=no_sord !for Strahler ordering
       else if(ORDER_SYSTEM(1:5).eq.'horsf')then
          nindex = no_hord !for Horsfield ordering
       endif

       ! get the inlet element 
       ! if two inlets, get the one with a higher order
       start_elem = umbilical_inlets(1)
       if(count(umbilical_inlets.NE.0).GT.1)then
           do indx=1,2
              inlet_ord(indx) = elem_ordrs(nindex,umbilical_inlets(indx))
           enddo
           if(inlet_ord(1).LT.inlet_ord(2))then
              start_elem = umbilical_inlets(2)
           endif
       endif

       start_rad = elem_field(ne_radius,start_elem)
       n_max_ord = elem_ordrs(nindex,start_elem)

       if(diagnostics_level.GT.1)then
          print *, "order system = ",order_system
	  print *, "s_ratio=",s_ratio
          print *, "calculated radii:"
          print *, "start radius start_rad = ",start_rad
       endif

       do ne=minelem_no_radius,remaining_elems
          radius=10.0_dp**(log10(s_ratio)*dble(elem_ordrs(nindex,ne)-n_max_ord)&
           +log10(start_rad))
          elem_field(ne_radius,ne)=radius  
          elem_field(ne_radius_in,ne)=radius
          elem_field(ne_radius_out,ne)=radius
          if(diagnostics_level.GT.1)then
	     print *,"element order for element",ne,"=",elem_ordrs(nindex,ne)
     	     print *,"radius for element",ne,"=",elem_field(ne_radius,ne)
     	     print *,"radius in for element",ne,"=",elem_field(ne_radius_in,ne)
     	     print *,"radius out for element",ne,"=",elem_field(ne_radius_out,ne)
          endif
       enddo

    endif !(minelem_no_radius.LT.num_elems)


    !calculate element volume
    do ne=1,num_elems
       elem_field(ne_vol,ne) = PI * elem_field(ne_radius,ne)**2 * &
            elem_field(ne_length,ne)
    enddo
    
    call enter_exit(sub_name,2)

  END subroutine define_rad_from_file
!
!##################################################################################
!

  subroutine define_rad_from_geom(ORDER_SYSTEM, CONTROL_PARAM, START_FROM, START_RAD, group_type_in, group_option_in)
  !*Description:* Defines vessel radius based on their geometric structure
    use arrays,only: dp,num_elems,elem_field,elem_ordrs,maxgen,elem_cnct,num_arterial_elems            
    use indices
    use diagnostics, only: enter_exit,get_diagnostics_level
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_RAD_FROM_GEOM" :: DEFINE_RAD_FROM_GEOM
    
   character(LEN=100), intent(in) :: ORDER_SYSTEM,START_FROM
   character(LEN=100), optional :: group_type_in, group_option_in
   real(dp), intent(in) :: CONTROL_PARAM,START_RAD
   !Input options ORDER_SYSTEM=STRAHLER (CONTROL_PARAM=RDS), HORSFIELD (CONTROL_PARAM=RDH)

   !Local variables
   character(LEN=100) :: group_type, group_options
   integer :: ne_min,ne_max,nindex,ne,n_max_ord,n,ne_start,&
      inlet_count
   real(dp) :: radius
   character(len=60) :: sub_name
   integer :: diagnostics_level

   sub_name = 'define_rad_from_geom'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)
    !define list of elements you are going to operate on
    if(present(group_type_in))then
      group_type = group_type_in
    else!default to all
      group_type='all'
    endif
    if(group_type.eq.'all')then
       ne_min=1
       ne_max=num_elems
    elseif(group_type.eq.'efield')then
    elseif(group_type.eq.'arterial')then
       ne_min=1
       ne_max=num_arterial_elems
    elseif(group_type.eq.'venous')then
       ne_min=num_arterial_elems + 1
       ne_start=ne_min
       ne_max=num_elems          
    elseif(group_type.eq.'list')then
      read (START_FROM,'(I10)') ne_min
      read (group_option_in,'(I10)') ne_max
    endif
    if(diagnostics_level.GT.1)then
    	  print *,"ne_min=", ne_min, "ne_max=",ne_max
    	endif
    !Define start element
    if(group_type.ne.'venous')then
       if(START_FROM.eq.'inlet')then
          inlet_count=0
      	  do ne=ne_min,ne_max
             if(elem_cnct(-1,0,ne).eq.0)then
                inlet_count=inlet_count+1
                ne_start=ne
             endif
      	  enddo
      	  if(inlet_count.gt.1)then
             WRITE(*,*) ' More than one inlet in this group, using last found, ne = ',ne_start
          endif
       else!element number defined
          read (START_FROM,'(I10)') ne_start
       endif
    endif
    if(diagnostics_level.GT.1)then
       print *, "ne_start=",ne_start
    endif

    !Strahler and Horsfield ordering system
    if(ORDER_SYSTEM(1:5).EQ.'strah')THEN
      nindex = no_sord !for Strahler ordering
    else if(ORDER_SYSTEM(1:5).eq.'horsf')then
      nindex = no_hord !for Horsfield ordering
    endif

    if(diagnostics_level.GT.1)then
       print *, "ordering system= ",ORDER_SYSTEM(1:5)
    endif

    ne=ne_start
    n_max_ord=elem_ordrs(nindex,ne)
    elem_field(ne_radius,ne)=START_RAD

    if(diagnostics_level.GT.1)then
       print *, "start radius START_RAD = ",START_RAD
    endif

    do ne=ne_min,ne_max
     radius=10.0_dp**(log10(CONTROL_PARAM)*dble(elem_ordrs(nindex,ne)-n_max_ord)&
        +log10(START_RAD))
     elem_field(ne_radius,ne)=radius  
     elem_field(ne_radius_in,ne)=radius
     elem_field(ne_radius_out,ne)=radius

     !element volume
     elem_field(ne_vol,ne) = PI * elem_field(ne_radius,ne)**2 * &
            elem_field(ne_length,ne)


     if(diagnostics_level.GT.1)then
	print *,"element order for element",ne,"=",elem_ordrs(nindex,ne)
     	print *,"radius for element",ne,"=",elem_field(ne_radius,ne)
     	print *,"radius in for element",ne,"=",elem_field(ne_radius_in,ne)
     	print *,"radius out for element",ne,"=",elem_field(ne_radius_out,ne)
        print *,"volume for element",ne,"=",elem_field(ne_vol,ne)
     endif
    enddo

    call enter_exit(sub_name,2)

  END subroutine define_rad_from_geom
!
!###########################################################################
!
  subroutine define_ven_rad_from_art(FILENAME,factor)
  !*Description:* Sets venous radii to corresponding arterial radii * factor
  !   (default factor is 2). 
  ! Prints arterial and corresponding venous radii to a file.  

    use arrays,only: dp,elem_field,num_elems,num_arterial_elems, elem_ordrs,&
                 art_ven_elem_map
    use indices,only: ne_length,ne_radius,no_sord,ne_vol 
    use other_consts, only: MAX_FILENAME_LEN
    use diagnostics, only: enter_exit,get_diagnostics_level
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_VEN_RAD_FROM_ART" :: DEFINE_VEN_RAD_FROM_ART

    character(len=MAX_FILENAME_LEN), intent(in) :: FILENAME
    real(dp), optional :: factor

    !     Local Variables
    real(dp), allocatable :: node_radius(:)
    integer :: ne,ne_art,ne_ven,art_sord,ven_sord
    real(dp) :: art_radius,ven_radius,rad_difference
    character(len=MAX_STRING_LEN) ::  venous_option
    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'define_ven_rad_from_art'
    call enter_exit(sub_name,1)
    
    !set venous radii to corresponding arterial radii * factor
    do ne_art=1,num_arterial_elems
       ne_ven = art_ven_elem_map(ne_art)
       if(ne_ven.GT.0)then
          art_radius = elem_field(ne_radius,ne_art)
          elem_field(ne_radius,ne_ven) = elem_field(ne_radius,ne_art) * factor
       endif
    enddo 

    !calculate element volume
    do ne=num_arterial_elems+1,num_elems
       elem_field(ne_vol,ne) = PI * elem_field(ne_radius,ne)**2 * &
            elem_field(ne_length,ne)
    enddo

    !print arterial and venous radii to a file 
    open(10, file=FILENAME, status="replace")
    write(10,*) 'arterial_elem,art_Strahler_order,arterial_elem_radius,venous_elem,&
                      ven_Strahler_order,venous_elem_radius,difference'
    !arterial elements and their matching venous elements
    do ne_art=1,num_arterial_elems
       ne_ven = art_ven_elem_map(ne_art)
       art_radius = elem_field(ne_radius,ne_art)
       art_sord = elem_ordrs(no_sord,ne_art)
       if(ne_ven.GT.0)then
          ven_radius = elem_field(ne_radius,ne_ven)
          ven_sord = elem_ordrs(no_sord,ne_ven)
          rad_difference = art_radius - ven_radius
          write(10,*) ne_art,',',art_sord, ',',art_radius,',',ne_ven,',',&
                       ven_sord,',',ven_radius,',',rad_difference
       else
          write(10,*) ne_art, ',',art_sord, ',', art_radius,',', 0, ',', 0,',', 0.0_dp, ',', 0.0_dp
       endif
    enddo
    close(10)
  
    call enter_exit(sub_name,2)

  END subroutine define_ven_rad_from_art
!
!##################################################################################
!
  subroutine element_connectivity_1d()
  !*Description:* Calculates element connectivity in 1D and stores in elem_cnct
    use arrays,only: elem_cnct,elem_nodes,elems_at_node,num_elems,num_nodes,elem_cnct_no_anast,&
                     anastomosis_elem
    use diagnostics, only: enter_exit,get_diagnostics_level
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ELEMENT_CONNECTIVITY_1D" :: ELEMENT_CONNECTIVITY_1D
  
    !     Local Variables
    integer :: ne,ne2,nn,noelem,np,np2,np1,counter,orphan_counter,np2_1
    integer,parameter :: NNT=2
    character(len=60) :: sub_name
    integer :: orphan_nodes(num_nodes)
    integer :: diagnostics_level

    sub_name = 'element_connectivity_1d'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    ! calculate elems_at_node array: stores the elements that nodes are in
    ! elems_at_node(node np,0)= total number of elements connected to this node
    ! elems_at_node(node np, index of each connected element starting at 1) = connected element
    elems_at_node = 0 !initialise

    DO ne=1,num_elems
       DO nn=1,2
          np=elem_nodes(nn,ne)
          elems_at_node(np,0)=elems_at_node(np,0)+1
          elems_at_node(np,elems_at_node(np,0))=ne ! local element that np is in
        ENDDO !nn
    ENDDO !noelem

    if(diagnostics_level.GT.1)then
    		DO nn=1,num_nodes
       		print *," "
       		print *,"node",nn
       		print *,"total number of elements connected",elems_at_node(nn,0)
       		DO ne=1,elems_at_node(nn,0)
          		print *,"element",elems_at_node(nn,ne)
       		ENDDO
    		ENDDO
    endif
    
    !check for nodes with 0 elements - exit if any are found
    orphan_counter = 0
    DO nn=1,num_nodes
		if(elems_at_node(nn,0).EQ.0)then
			orphan_counter = orphan_counter + 1
			orphan_nodes(orphan_counter) = nn
		endif
    ENDDO
    if(orphan_counter.GT.0)then
		print *, "found",orphan_counter,"node(s) not connected to any elements"
		do counter=1,orphan_counter
			print *,"node",orphan_nodes(counter),"is not connected to any elements"
		enddo
		call exit(0)
    endif
    
    elem_cnct=0 !initialise all elem_cnct
    elem_cnct_no_anast = 0 !initialise

    DO ne=1,num_elems
       !     ne_global=elems(noelem)
       IF(NNT == 2) THEN !1d
          np1=elem_nodes(1,ne) !first local node
          np2=elem_nodes(2,ne) !second local node

          DO noelem=1,elems_at_node(np2,0) !for each element connected to node np2
             ne2=elems_at_node(np2,noelem) !get the element number connected to node np2
             IF(ne2 /= ne)THEN !if element connected to node np2 is not the current element ne
                !check that the second node of the current element is the first node of ne2
                np2_1 = elem_nodes(1,ne2)
                if(np2.EQ.np2_1)then
                   elem_cnct(-1,0,ne2)=elem_cnct(-1,0,ne2)+1
                   elem_cnct(-1,elem_cnct(-1,0,ne2),ne2)=ne !previous element              
                   elem_cnct(1,0,ne)=elem_cnct(1,0,ne)+1
                   elem_cnct(1,elem_cnct(1,0,ne),ne)=ne2
                   if((ne2.NE.anastomosis_elem).AND.(ne.NE.anastomosis_elem))then
                      elem_cnct_no_anast(-1,0,ne2)=elem_cnct_no_anast(-1,0,ne2)+1
                      elem_cnct_no_anast(-1,elem_cnct_no_anast(-1,0,ne2),ne2)=ne !previous element              
                      elem_cnct_no_anast(1,0,ne)=elem_cnct_no_anast(1,0,ne)+1
                      elem_cnct_no_anast(1,elem_cnct_no_anast(1,0,ne),ne)=ne2
                   endif
                endif
             ENDIF !ne2
          ENDDO !noelem2


       ENDIF
    ENDDO

	! total count of upstream elements connected to element ne elem_cnct(-1,0,ne)
	! upstream elements elem_cnct(-1,counter,ne)
	! total count of downstream elements connected to element ne elem_cnct(1,0,ne)
	! downstream elements elem_cnct(1,counter,ne)
    if(diagnostics_level.GT.1)then
   		DO ne=1,num_elems
   	    		print *,""
   	    		print *,"element",ne 
       		IF(elem_cnct(-1,0,ne).gt.0)THEN
       	    		print *,"total number of upstream elements:",elem_cnct(-1,0,ne)
       			DO counter=1,elem_cnct(-1,0,ne)
          			print *,"upstream element",elem_cnct(-1,counter,ne)
       	    		ENDDO
       		ENDIF
       		IF(elem_cnct(1,0,ne).gt.0)THEN
       	    		print *,"total number of downstream elements:",elem_cnct(1,0,ne)
       			DO counter=1,elem_cnct(1,0,ne)
          			print *,"downstream element",elem_cnct(1,counter,ne)
       	    		ENDDO
       		ENDIF
    		ENDDO
                print *,"element connectivity without the anastomosis:"
  		DO ne=1,num_elems
   	    		print *,""
   	    		print *,"element",ne 
       		IF(elem_cnct_no_anast(-1,0,ne).gt.0)THEN
       	    		print *,"total number of upstream elements:",elem_cnct_no_anast(-1,0,ne)
       			DO counter=1,elem_cnct_no_anast(-1,0,ne)
          			print *,"upstream element",elem_cnct_no_anast(-1,counter,ne)
       	    		ENDDO
       		ENDIF
       		IF(elem_cnct_no_anast(1,0,ne).gt.0)THEN
       	    		print *,"total number of downstream elements:",elem_cnct_no_anast(1,0,ne)
       			DO counter=1,elem_cnct_no_anast(1,0,ne)
          			print *,"downstream element",elem_cnct_no_anast(1,counter,ne)
       	    		ENDDO
       		ENDIF
    		ENDDO
    endif

    call enter_exit(sub_name,2)

  END subroutine element_connectivity_1d
!
!###################################################################################
!
  subroutine evaluate_ordering()
  !*Description:* calculates generations, Horsfield orders, Strahler orders for a given tree
   
    use arrays,only: elem_cnct_no_anast,elem_nodes,elem_ordrs,elem_symmetry,&
         elems_at_node,num_elems,num_nodes,maxgen,anastomosis_elem
    use diagnostics, only: enter_exit,get_diagnostics_level
    implicit none
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_ORDERING" :: EVALUATE_ORDERING

    integer :: INLETS,ne,ne0,ne2,noelem2,np,np2, &
         num_attach,n_children,n_generation, &
         n_horsfield,OUTLETS,STRAHLER,STRAHLER_ADD,temp1, counter
    LOGICAL :: DISCONNECT,DUPLICATE
    character(len=60) :: sub_name
    integer :: diagnostics_level
    
    sub_name = 'evaluate_ordering'    
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    !Calculate generations, Horsfield orders, Strahler orders
    !.....Calculate branch generations

    elem_ordrs = 0
    maxgen=1
    DO ne=1,num_elems
       ne0=elem_cnct_no_anast(-1,1,ne) !parent
       IF(ne0.NE.0)THEN
          n_generation=elem_ordrs(1,ne0) !parent generation
          IF(elem_cnct_no_anast(1,0,ne0).EQ.1)THEN !single daughter
             elem_ordrs(1,ne)=n_generation + (elem_symmetry(ne)-1)
          ELSE IF(elem_cnct_no_anast(1,0,ne0).GE.2)THEN
             elem_ordrs(1,ne)=n_generation+1
          ENDIF
       ELSE
          elem_ordrs(1,ne)=1 !generation 1
       ENDIF
       maxgen=max(maxgen,elem_ordrs(1,ne))

    ENDDO !noelem

    if(diagnostics_level.GT.1)then
	print *,"branch generations - maxgen",maxgen
    endif

    !.....Calculate the branch orders
    DO ne=num_elems,1,-1
       n_horsfield=MAX(elem_ordrs(2,ne),1)
       n_children=elem_cnct_no_anast(1,0,ne) !number of child branches
       IF(n_children.EQ.1)THEN
          IF(elem_ordrs(1,elem_cnct_no_anast(1,1,ne)).EQ.0)  n_children=0
       ENDIF
       STRAHLER=0
       STRAHLER_ADD=1
       IF(n_children.GE.2)THEN !branch has two or more daughters
          STRAHLER=elem_ordrs(3,elem_cnct_no_anast(1,1,ne)) !first daughter
          DO noelem2=1,n_children !for all daughters
             ne2=elem_cnct_no_anast(1,noelem2,ne) !global element # of daughter
             temp1=elem_ordrs(2,ne2) !Horsfield order of daughter
             IF(temp1.GT.n_horsfield) n_horsfield=temp1
             IF(elem_ordrs(3,ne2).LT.STRAHLER)THEN
                STRAHLER_ADD=0
             ELSE IF(elem_ordrs(3,ne2).GT.STRAHLER)THEN
                STRAHLER_ADD=0
                STRAHLER=elem_ordrs(3,ne2) !highest daughter
             ENDIF
          ENDDO !noelem2 (ne2)
          n_horsfield=n_horsfield+1 !Horsfield ordering
       ELSE IF(n_children.EQ.1)THEN
          ne2=elem_cnct_no_anast(1,1,ne) !local element # of daughter
          n_horsfield=elem_ordrs(2,ne2)+(elem_symmetry(ne)-1)
          STRAHLER_ADD=elem_ordrs(3,ne2)+(elem_symmetry(ne)-1)
       ENDIF !elem_cnct_no_anast
       elem_ordrs(2,ne)=n_horsfield !store the Horsfield order
       elem_ordrs(3,ne)=STRAHLER+STRAHLER_ADD !Strahler order
    ENDDO !noelem

    !if anastomosis element exists set its orders to 0
    if(anastomosis_elem.GT.0)then
	elem_ordrs(1,anastomosis_elem)=0
	elem_ordrs(2,anastomosis_elem)=0
	elem_ordrs(3,anastomosis_elem)=0
    endif


    !       Check for disconnected nodes and number of inlets and outlets
    DUPLICATE=.FALSE.
    DO ne=1,num_elems
       np=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       IF(np.EQ.np2)THEN
          DUPLICATE=.TRUE.
       ENDIF
    ENDDO

    DISCONNECT=.FALSE.
    INLETS=0
    OUTLETS=0
    DO np=1,num_nodes
       num_attach=elems_at_node(np,0)
       IF(num_attach.EQ.0)THEN
          DISCONNECT=.TRUE.
       ELSEIF(num_attach.EQ.1)THEN
          ne=elems_at_node(np,1)
          IF(elem_cnct_no_anast(1,0,ne).EQ.0) OUTLETS=OUTLETS+1
         IF(elem_cnct_no_anast(-1,0,ne).EQ.0) INLETS=INLETS+1
       ELSEIF(num_attach.GT.3)THEN
          WRITE(*,*) ' Node ',np,' attached to',num_attach,' elements'
       ENDIF
    ENDDO
  
    do ne=1,num_elems
       if(elem_ordrs(3,ne).EQ.0)then
          if((anastomosis_elem.EQ.0).OR.((anastomosis_elem.GT.0).AND.(ne.NE.anastomosis_elem)))then
             print *, "Warning: element", ne,"was assigned Strahler order 0"
          endif
       endif
    enddo

    if(diagnostics_level.GT.1)then 
       do ne=1,num_elems
          do counter=1,3		
             print *, "elem_ordrs(",counter,",",ne,")=",elem_ordrs(counter,ne)
    	  enddo
       enddo	  
    endif

    call enter_exit(sub_name,2)

  end subroutine evaluate_ordering
!
!###################################################################################
!
  subroutine reallocate_node_elem_arrays(num_elems_new,num_nodes_new)
  !*Description:* Reallocates the size of arrays when modifying geometries
  
    use arrays,only: dp,elems,elem_cnct,elem_direction,elem_field,&
         elem_ordrs,elem_nodes,&
         elem_symmetry,elem_units_below,elems_at_node,&
         nodes,node_field,node_xyz,num_elems,num_nodes
    use indices
    use diagnostics, only: enter_exit
    implicit none

!!! Parameters
    integer,intent(in) :: num_elems_new,num_nodes_new

!!! Local variables
    integer,allocatable :: nodelem_temp(:),enodes_temp(:,:),enodes_temp2(:,:,:)
    real(dp),allocatable :: xyz_temp(:,:),rnodes_temp(:,:)
    character(len=60) :: sub_name

    sub_name = 'reallocate_node_elem_arrays'
    call enter_exit(sub_name,1)

    allocate(nodelem_temp(num_nodes))
    nodelem_temp = nodes ! copy to temporary array
    deallocate(nodes) !deallocate initially allocated memory
    allocate(nodes(num_nodes_new))
    nodes(1:num_nodes)=nodelem_temp(1:num_nodes)
    deallocate(nodelem_temp) !deallocate the temporary array

    allocate(xyz_temp(3,num_nodes))
    xyz_temp=node_xyz
    deallocate(node_xyz)
    allocate(node_xyz(3,num_nodes_new))
    node_xyz(1:3,1:num_nodes)=xyz_temp(1:3,1:num_nodes)

    allocate(nodelem_temp(num_elems))
    nodelem_temp = elems ! copy to temporary array
    deallocate(elems) !deallocate initially allocated memory
    allocate(elems(num_elems_new))
    elems(1:num_elems)=nodelem_temp(1:num_elems)
    deallocate(nodelem_temp) !deallocate the temporary array

    allocate(enodes_temp(2,num_elems))
    enodes_temp=elem_nodes
    deallocate(elem_nodes)
    allocate(elem_nodes(2,num_elems_new))
    elem_nodes(1:2,1:num_elems)=enodes_temp(1:2,1:num_elems)
    deallocate(enodes_temp)

    allocate(rnodes_temp(num_ne,num_elems))
    rnodes_temp=elem_field
    deallocate(elem_field)
    allocate(elem_field(num_ne,num_elems_new))
    elem_field(1:num_ne,1:num_elems)=rnodes_temp(1:num_ne,1:num_elems)
    deallocate(rnodes_temp)
    elem_field(1:num_ne,num_elems+1:num_elems_new) = 0.0_dp

    allocate(rnodes_temp(3,num_elems))
    rnodes_temp=elem_direction
    deallocate(elem_direction)
    allocate(elem_direction(3,num_elems_new))
    elem_direction(1:3,1:num_elems)=rnodes_temp(1:3,1:num_elems)
    deallocate(rnodes_temp)
    elem_direction(1:3,num_elems+1:num_elems_new) = 0.0_dp

    allocate(rnodes_temp(num_nj,num_nodes))
    rnodes_temp=node_field
    deallocate(node_field)
    allocate(node_field(num_nj,num_nodes_new))
    node_field(1:num_nj,1:num_nodes)=rnodes_temp(1:num_nj,1:num_nodes)
    deallocate(rnodes_temp)
    node_field(1:num_nj,num_nodes+1:num_nodes_new)=0.0_dp

    allocate(nodelem_temp(num_elems))
    nodelem_temp = elem_symmetry ! copy to temporary array
    deallocate(elem_symmetry) !deallocate initially allocated memory
    allocate(elem_symmetry(num_elems_new))
    elem_symmetry(1:num_elems)=nodelem_temp(1:num_elems)
    deallocate(nodelem_temp) !deallocate the temporary array
    elem_symmetry(num_elems+1:num_elems_new)=1

    allocate(enodes_temp2(-1:1,0:2,0:num_elems))
    enodes_temp2=elem_cnct
    deallocate(elem_cnct)
    allocate(elem_cnct(-1:1,0:2,0:num_elems_new))
    elem_cnct(-1:1,0:2,0:num_elems)=enodes_temp2(-1:1,0:2,0:num_elems)
    deallocate(enodes_temp2)
    elem_cnct(-1:1,0:2,num_elems+1:num_elems_new) = 0

    allocate(enodes_temp(num_ord,num_elems))
    enodes_temp=elem_ordrs
    deallocate(elem_ordrs)
    allocate(elem_ordrs(num_ord,num_elems_new))
    elem_ordrs(1:num_ord,1:num_elems)=enodes_temp(1:num_ord,1:num_elems)
    deallocate(enodes_temp)
    elem_ordrs(1:num_ord,num_elems+1:num_elems_new) = 0

    allocate(nodelem_temp(num_elems))
    nodelem_temp=elem_units_below
    deallocate(elem_units_below)
    allocate(elem_units_below(num_elems_new))
    elem_units_below(1:num_elems)=nodelem_temp(1:num_elems)
    deallocate(nodelem_temp)
    elem_units_below(num_elems+1:num_elems_new)=0

    allocate(enodes_temp(num_nodes,0:3))
    enodes_temp=elems_at_node
    deallocate(elems_at_node)
    allocate(elems_at_node(num_nodes_new,0:3))
    elems_at_node(1:num_nodes,0:3)=enodes_temp(1:num_nodes,0:3)
    deallocate(enodes_temp)
    elems_at_node(num_nodes+1:num_nodes_new,0:3)=0

    call enter_exit(sub_name,2)

  end subroutine reallocate_node_elem_arrays

!
!###################################################################################
!

subroutine update_1d_elem_field(ne_field,elem_number,value)

    use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
    use arrays,only: dp,elem_field,num_elems
    use indices,only: num_ne
    use diagnostics, only: enter_exit, get_diagnostics_level
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_UPDATE_1D_ELEM_FIELD" :: UPDATE_1D_ELEM_FIELD

!!! Parameters
    integer, intent(in) :: ne_field
    integer, intent(in) :: elem_number
    real(dp), intent(in) :: value

    character(len=60) :: sub_name
    integer :: diagnostics_level

    sub_name = 'update_1d_elem_field'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    if(ne_field.gt.num_ne)then
      print *, 'WARNING: trying to allocate to elem field greater than defined for this problem',ne_field,num_ne
    elseif(elem_number.gt.num_elems)then
      print *, 'WARNING: trying to allocate to element greater than defined for this problem',elem_number,num_elems
    else
      elem_field(ne_field,elem_number) = value
    endif
    call enter_exit(sub_name,2)
  end subroutine update_1d_elem_field
 !
 !##############################################################
 !

  subroutine get_final_real(string,rtemp)
    use arrays,only: dp
    implicit none
    character, intent(in) :: string*(132)
    integer :: ibeg,iend
    real(dp), intent(out) :: rtemp
    real(dp) :: rsign
    character :: sub_string*(40)

    iend=len(string) !get the length of the string
    ibeg=index(string,":")+1 !get location of real in string
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond :
    iend=len(sub_string) !get the length of the sub-string
    if(sub_string(1:1).eq.'-')then !check whether negative
       rsign=-1.0_dp
       ibeg=2
    else
       rsign=1.0_dp
       ibeg=1
    endif
    read (sub_string(ibeg:iend), * ) rtemp !get real value
    rtemp=rtemp*rsign !apply sign to number

  end subroutine get_final_real

!
!###################################################################################
!
  subroutine get_final_string(string,rtemp)
    use arrays,only: dp
    implicit none
    character, intent(in) :: string*(132)
    integer :: ibeg,iend
    real(dp), intent(out) :: rtemp
    real(dp) :: rsign
    character :: sub_string*(40)

    iend=len(string) !get the length of the string
    ibeg=index(string,":")+1 !get location of real in string
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond :
    iend=len(sub_string) !get the length of the sub-string
    if(sub_string(1:1).eq.'-')then !check whether negative
       rsign=-1.0_dp
       ibeg=2
    else
       rsign=1.0_dp
       ibeg=1
    endif
    read (sub_string(ibeg:iend), '(D25.17)' ) rtemp !get real value
    rtemp=rtemp*rsign !apply sign to number

  end subroutine get_final_string

!
!###################################################################################
!
  subroutine get_local_node(np_global,np_local)
    use arrays,only: nodes,num_nodes
    implicit none

    integer,intent(in) :: np_global
    integer,intent(out) :: np_local

    integer :: np
    logical :: found

    np=1   
    found=.false.
    do while ((.not.found).AND.(np.le.num_nodes))
       if(nodes(np).eq.np_global)then
          found=.true.
       else
          np=np+1
       endif
    enddo

    if(.not.found)then
       np = 0          
       write(*,'('' Global node '',I6,'' not in node list'')') np_global
    endif

    np_local = np
    return

  end subroutine get_local_node

!
!###################################################################################
!
  subroutine get_final_integer(string,num)
    implicit none
    character,intent(in) :: string*(132)
    integer,intent(out) :: num
    integer :: ibeg,iend,nsign,ntemp
    character :: sub_string*(40)

    iend=len(string) !get the length of the string
    ibeg=index(string,":")+1 !get location of integer in string, follows ":"
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond ":"
    iend=len(sub_string) !length of the sub-string
    if(sub_string(1:1).eq.'-')then !check for negative sign
       nsign=-1
       ibeg=2
    else
       nsign=1
       ibeg=1
    endif
    read (sub_string(ibeg:iend), '(i10)' ) ntemp !get integer values
    ntemp=ntemp*nsign !apply sign to number

    num=ntemp !return the integer value

  end subroutine get_final_integer

!
! ##########################################################################      
!

  function inlist(item,ilist)
!!! dummy arguments
    integer :: item,ilist(:)
! local variables
    integer :: n
    logical :: inlist

    inlist = .false.
    do n=1,size(ilist)
       if(item == ilist(n)) inlist = .true.
    enddo

  end function inlist
!
!###########################################################################################
!
end module geometry

