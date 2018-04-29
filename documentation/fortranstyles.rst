
========================
Fortran style guidelines
========================

Here are some examples of fortran programming style to follow 

Preamble
========
::

    module module name
    !*Description:* description of the module 
    !
    !*Contributor(s): list the contributor names
    !
      use anyothermodulesused
      implicit none
  
      !Module parameters
  
      !Module types

      !Module variables

      !Interfaces
      private list1
      public list2
      public list3

    contains
    !
    !##############################################################################
    !

Indenting
=========
Please indent using two spaces.

Subroutines
===========
Include subroutines in alphabetical order.

Separate subroutines as follows
::

    !
    !##############################################################################
    !

Begin each subroutine with a brief description of its functionality.
This will be automatically included in the documentation. The description must follow
the subroutine opening statement: subroutine subname(things,you,pass)
::

    subroutine subname(things,you,pass)
    !*Description:* What the subroutine does
    
      use module1,only: stuff,from,module1 
      use module2
      use diagnostics, only: enter_exit
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SUBNAME: SUBNAME
    
      integer, intent(in):: things
      integer, intent(out) :: you
      real(dp), intent(inout)  :: pass
      !local variables
      integer :: other,stuff

      character(len=60) :: sub_name

      sub_name = 'subname'
      call enter_exit(sub_name,1)

      !BODY OF CODE

      call enter_exit(sub_name,2)
    end subroutine subname

Note that if a subroutine has many arguments on multiple lines, the description of the 
subroutine won't be included in the generated documentation. In that case add the subroutine 
description to the module description.

Functions
=========
Include functions in alphabetical order below subroutines.
