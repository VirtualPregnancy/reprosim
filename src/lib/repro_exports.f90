!
!*Description:* This module handles all export functions
!
module repro_exports
  implicit none

  private
  public export_1d_elem_geometry,export_node_geometry,export_node_field,&
       export_terminal_perfusion,&
       export_1d_elem_field

contains
!!!################################################################

  subroutine export_1d_elem_field(ne_field, EXELEMFILE, group_name, field_name )
 
    use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
    use arrays,only: elem_field,num_elems
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_1D_ELEM_FIELD" :: EXPORT_1D_ELEM_FIELD
  
!!! Parameters
    integer, intent(in) :: ne_field
    character(len=MAX_FILENAME_LEN), intent(in) :: EXELEMFILE
    character(len=MAX_STRING_LEN), intent(in) :: field_name
    character(len=MAX_STRING_LEN), intent(in) :: group_name

!!! Local Variables
    integer :: len_end,ne
    logical :: CHANGED

    open(10, file=EXELEMFILE, status='replace')

    len_end=len_trim(group_name)
    !**     write the group name
    write(10,'( '' Group name: '',A)') group_name(:len_end)
    !**         write the elements
    write(10,'( '' Shape.  Dimension=1'' )')
    CHANGED=.TRUE. !initialise to force output of element information
    len_end=len_trim(field_name)
    do ne=1,num_elems
       if(ne>1) THEN
          CHANGED=.FALSE.
       endif
       if(CHANGED)THEN
          write(10,'( '' #Scale factor sets=0'' )')
          write(10,'( '' #Nodes= 0'' )')
          write(10,'( '' #Fields= 1'' )')
          write(10,'( '' 1)'',A,'', field, rectangular cartesian, #Components=1'')')&
               field_name(:len_end)
          write(10,'( ''  '',A,''.  l.Lagrange, no modify, grid based.'')') &
               field_name(:len_end)
          write(10,'( ''  #xi1=1'')')
       endif

       write(10,'(1X,''Element: '',I12,'' 0 0'' )') ne
       write(10,'(3X,''Values:'' )')
       write(10,'(4X,2(1X,E12.5))') elem_field(ne_field,ne),elem_field(ne_field,ne)
    enddo !no_nelist (ne)
    close(10)

  end subroutine export_1d_elem_field

!!!############################################################################

  subroutine export_1d_elem_geometry(EXELEMFILE, name)

    use arrays,only: elem_nodes,num_elems
    use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
    use diagnostics, only: enter_exit

    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_1D_ELEM_GEOMETRY" :: EXPORT_1D_ELEM_GEOMETRY

!!! Parameters
    character(len=MAX_FILENAME_LEN), intent(in) :: EXELEMFILE
    character(len=MAX_STRING_LEN), intent(in) :: name

!!! Local Variables
    integer :: len_end,ne,nj,nn
    character(len=1) :: char1
    logical :: CHANGED
    character(len=60) :: sub_name

    sub_name = 'export_1d_elem_geometry'
    call enter_exit(sub_name,1)

    open(10, file=EXELEMFILE, status='replace')
    len_end=len_trim(name)
    !**     write the group name
    write(10,'( '' Group name: '',A)') name(:len_end)
    !**         write the elements
    write(10,'( '' Shape.  Dimension=1'' )')
    CHANGED=.TRUE. !initialise to force output of element information
    do ne=1,num_elems
       if(ne>1) THEN
          CHANGED=.FALSE.
       endif
       if(CHANGED)THEN
          write(10,'( '' #Scale factor sets=1'' )')
          write(10,'( ''   l.Lagrange, #Scale factors= 2'' )')
          write(10,'( '' #Nodes= 2'' )')
          write(10,'( '' #Fields= 1'' )')
          write(10,'( '' 1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
          do nj=1,3
             if(nj==1) char1='x'; if(nj==2) char1='y'; if(nj==3) char1='z';
             write(10,'(''  '',A2,''.  l.Lagrange, no modify, standard node based.'')') char1
             write(10,'( ''     #Nodes= 2'')')
             do nn=1,2
                write(10,'(''      '',I1,''.  #Values=1'')') nn
                write(10,'(''       Value indices:      1 '')')
                write(10,'(''       Scale factor indices:'',I4)') nn
             enddo !nn
          enddo !nj
       endif
       write(10,'(1X,''Element: '',I12,'' 0 0'' )') ne
       !**               write the nodes
       write(10,'(3X,''Nodes:'' )')
       write(10,'(4X,2(1X,I12))') elem_nodes(1,ne),elem_nodes(2,ne)
       !**                 write the scale factors
       write(10,'(3X,''Scale factors:'' )')
       write(10,'(4X,2(1X,E12.5))') 1.d0,1.d0
    enddo !no_nelist (ne)
    close(10)

    call enter_exit(sub_name,2)
  end subroutine export_1d_elem_geometry


!!!##########################################################################

  subroutine export_node_geometry(EXNODEFILE, name)
  
    use arrays,only: node_xyz,num_nodes
    use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_NODE_GEOMETRY" :: EXPORT_NODE_GEOMETRY

!!! Parameters
    character(len=MAX_FILENAME_LEN),intent(in) :: EXNODEFILE
    character(len=MAX_STRING_LEN),intent(in) :: name

!!! Local Variables
    integer :: len_end,nj,np,np_last,VALUE_INDEX
    logical :: FIRST_NODE

    len_end=len_trim(name)
    if(num_nodes.GT.0) THEN
       open(10, file=EXNODEFILE, status='replace')
       !**     write the group name
       write(10,'( '' Group name: '',A)') name(:len_end)
       FIRST_NODE=.TRUE.
       np_last=1
       !*** Exporting Geometry
       do np=1,num_nodes
          if(np.gt.1) np_last = np
          !*** Write the field information
          VALUE_INDEX=1
          if(FIRST_NODE)THEN
             write(10,'( '' #Fields=1'' )')
             write(10,'('' 1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
             do nj=1,3
                if(nj.eq.1) write(10,'(2X,''x.  '')',advance="no")
                if(nj.eq.2) write(10,'(2X,''y.  '')',advance="no")
                if(nj.eq.3) write(10,'(2X,''z.  '')',advance="no")
                write(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="no") 1,0
                write(10,'()')
             enddo
          endif !FIRST_NODE
          !***      write the node
          write(10,'(1X,''Node: '',I12)') np
          do nj=1,3
             write(10,'(2X,4(1X,F12.6))') (node_xyz(nj,np))
          enddo !njj2
          FIRST_NODE=.FALSE.
          np_last=np
       enddo !nolist (np)
    endif !num_nodes
    close(10)

  end subroutine export_node_geometry

!!!########################################################################

  subroutine export_terminal_perfusion(EXNODEFILE, name)
  
    use arrays,only: elem_nodes,&
         node_xyz,num_units,units,unit_field
    use indices
    use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_TERMINAL_PERFUSION" :: EXPORT_TERMINAL_PERFUSION

!!! Parameters
    character(len=MAX_FILENAME_LEN),intent(in) :: EXNODEFILE
    character(len=MAX_STRING_LEN),intent(in) :: name

!!! Local Variables
    integer :: len_end,ne,nj,NOLIST,np,np_last,VALUE_INDEX
    logical :: FIRST_NODE

    len_end=len_trim(name)
    if(num_units.GT.0) THEN
       open(10, file=EXNODEFILE, status='replace')
       !**     write the group name
       write(10,'( '' Group name: '',A)') name(:len_end)
       FIRST_NODE=.TRUE.
       np_last=1
       !*** Exporting Terminal Solution
       do nolist=1,num_units
          if(nolist.GT.1) np_last = np
          ne=units(nolist)
          np=elem_nodes(2,ne)
          !*** Write the field information
          VALUE_INDEX=1
          if(FIRST_NODE)THEN
             write(10,'( '' #Fields=3'' )')
             write(10,'('' 1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
             do nj=1,3
                if(nj.eq.1) write(10,'(2X,''x.  '')',advance="no")
                if(nj.eq.2) write(10,'(2X,''y.  '')',advance="no")
                if(nj.eq.3) write(10,'(2X,''z.  '')',advance="no")
                write(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
                VALUE_INDEX=VALUE_INDEX+1
             enddo
             !perfusion
             write(10,'('' 2) flow, field, rectangular cartesian, #Components=1'')')
             write(10,'(2X,''1.  '')',advance="no")
             write(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
             !Pressure
             VALUE_INDEX=VALUE_INDEX+1
             write(10,'('' 3) pressure, field, rectangular cartesian, #Components=1'')')
             write(10,'(2X,''1.  '')',advance="no")
             write(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") VALUE_INDEX,0
          endif !FIRST_NODE
          !***      write the node
          write(10,'(1X,''Node: '',I12)') np
          do nj=1,3
             write(10,'(2X,4(1X,F12.6))') (node_xyz(nj,np))      !Coordinates
          enddo !njj2
           write(10,'(2X,4(1X,F20.6))') (unit_field(nu_perf,NOLIST)) !flow
           write(10,'(2X,4(1X,F12.6))') (unit_field(nu_blood_press,NOLIST)) !pressure
          FIRST_NODE=.FALSE.
          np_last=np
       enddo !nolist (np)
    endif !num_nodes
    close(10)

  end subroutine export_terminal_perfusion




!!! #################################################################

  subroutine export_node_field(nj_field, EXNODEFIELD, name, field_name)
  
    use arrays,only: node_field,num_nodes
    use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
    implicit none
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EXPORT_NODE_FIELD" :: EXPORT_NODE_FIELD

!!! Parameters
    integer,intent(in) :: nj_field
    character(len=MAX_FILENAME_LEN),intent(in) :: EXNODEFIELD
    character(len=MAX_STRING_LEN),intent(in) :: field_name
    character(len=MAX_STRING_LEN),intent(in) :: name

!!! Local Variables
    integer :: len_end,np
    logical :: FIRST_NODE

    open(10, file=EXNODEFIELD, status='replace')
    !**     write the group name
    len_end=len_trim(name)
    write(10,'( '' Group name: '',A)') name(:len_end)
    len_end=len_trim(field_name)
    FIRST_NODE=.TRUE.
    !*** the field as specified by user
    do np=1,num_nodes
       !*** Write the field information
       if(FIRST_NODE)THEN
          write(10,'( '' #Fields=1'' )')
          write(10,'('' 1) '',A,'', field, rectangular cartesian, #Components=1'')') &
               field_name(:len_end)
          write(10,'(2X,''1.  '')',advance="no")
          write(10,'(''Value index='',I1,'', #Derivatives='',I1)',advance="yes") 1,0
       endif !FIRST_NODE
       !***      write the node
       write(10,'(1X,''Node: '',I12)') np
       write(10,'(2X,2(1X,F12.6))') (node_field(nj_field,np))
       FIRST_NODE=.FALSE.
    enddo !num_nodes
    close(10)

  end subroutine export_node_field

!!! ###########################################################

end module repro_exports
