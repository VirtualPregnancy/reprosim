module math_utilities
!*Description:* This module contains solvers required for blood vessel problems
!
  implicit none
  private
  public ax_cr,diagonal_pointer_cr,ilu_cr,lus_cr,mult_givens,rearrange_cr

contains
!
!###########################################################################
!
  subroutine ax_cr ( n, nz_num, ia, ja, a, x, w )
  !*Description:* Computes A*x for a matrix stored in sparse compressed row form
    use arrays, only: dp

    integer :: n !the order of the system
    integer :: nz_num
    integer :: ia(n+1) !row indices
    integer :: ja(nz_num) !column indices
    real(dp) :: a(nz_num) !Matrix values
    real(dp) :: x(n) !Vector to be multiplied by A
    real(dp) :: w(n) !Value of A*x

    integer :: i
    integer :: k1
    integer :: k2

    w(1:n) = 0.0D+00

    do i = 1, n
       k1 = ia(i)
       k2 = ia(i+1) - 1
       w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
    end do

    return
  end subroutine ax_cr
!
!##############################################################################
!
  subroutine ilu_cr ( n, nz_num, l_size, ia, ja, a, ua, l)
  ! *Description:* computes the incomplete LU factorization of a matrix. For a matrix
  ! stored in compressed row format.
  ! Input, integer :: UA(N), the index of the diagonal element of each row.
  ! Output, real(dp) :: L(NZ_NUM), the ILU factorization of A.
    use arrays, only: dp
    integer :: n
    integer :: nz_num
    integer :: l_size
    integer :: ia(n+1)
    integer :: ja(nz_num)
    real(dp) :: a(nz_num)
    integer:: ua(n)
    real(dp) :: l(l_size)

    integer :: i
    integer, allocatable :: iw(:)
    integer :: j
    integer :: jj
    integer :: jrow
    integer :: jw
    integer :: k
    real(dp) :: tl
    !  Copy A.
    l(1:nz_num) = a(1:nz_num)
    if(allocated(iw)) deallocate(iw)
    allocate(iw(n))
    do i = 1, n ! for each row, up to max number of rows
       !  IW points to the nonzero entries in row I.
       iw(1:n) = -1
       do k = ia(i), ia(i+1) - 1 !for each
          iw(ja(k)) = k
       end do
       do j = ia(i), ia(i+1) - 1
          jrow = ja(j)
          if ( i <= jrow ) then
             exit
          end if
          tl = l(j) * l(ua(jrow))
          l(j) = tl
          do jj = ua(jrow) + 1, ia(jrow+1) - 1
             jw = iw(ja(jj))
             if ( jw /= -1 ) then
                l(jw) = l(jw) - tl * l(jj)
             end if
          end do
       end do
       ua(i) = j
       if ( jrow /= i ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ILU_CR - Fatal error!'
          write ( *, '(a)' ) '  JROW ~= I'
          write ( *, '(a,i8)' ) '  JROW = ', jrow
          write ( *, '(a,i8)' ) '  I    = ', i
          stop
       end if
       if ( l(j) == 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ILU_CR - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
          write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
          stop
       end if
       l(j) = 1.0D+00 / l(j)
    end do

    l(ua(1:n)) = 1.0D+00 / l(ua(1:n))
    deallocate(iw)
    !return
  end subroutine ilu_cr
!
!##############################################################################
!
subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )
!*Description:* finds diagonal entries in a sparse compressed row matrix.
! The array UA can be used to locate the diagonal elements of the matrix.
! It is assumed that every row of the matrix includes a diagonal element,
! and that the elements of each row have been ascending sorted.
    use arrays, only: dp
    integer :: n
    integer :: nz_num
    integer :: ia(n+1)
    integer ::  ja(nz_num)
    integer :: ua(n)

    integer :: i
    integer :: k

    ua(1:n) = -1

    do i = 1, n
       do k = ia(i), ia(i+1) - 1
          if ( ja(k) == i ) then
             ua(i) = k
          end if
       end do
    end do
    return
  end subroutine diagonal_pointer_cr

  !*****************************************************************************80

  subroutine lus_cr ( n, nz_num,l_size, ia, ja, l, ua, r)
!*Description:* this subroutine applies the incomplete LU preconditioner.
! The linear system M * Z = R is solved for Z.  M is the incomplete
! LU preconditioner matrix, and R is a vector supplied by the user.
! So essentially, we're solving L * U * Z = R.
! Input, integer :: UA(N), the index of the diagonal element
! of each row.
! Input, real(dp) :: R(N), the right hand side.
! Output, real(dp) :: Z(N), the solution of the system M * Z = R.
    use arrays, only:dp

    integer :: n
    integer :: nz_num
    integer :: l_size
    integer :: ia(n+1)
    integer :: ja(nz_num)
    real(dp) :: l(l_size)
    integer :: ua(n)
    real(dp) :: r(n)

    integer :: i
    integer :: j
    real(dp), allocatable :: w(:)
    !real(dp) :: z(n)
    if(allocated(w)) deallocate(w)
    allocate(w(n))
    !  Copy R in.
    w(1:n) = r(1:n)

    !  Solve L * w = w where L is unit lower triangular.
    do i = 2, n
       do j = ia(i), ua(i) - 1
          w(i) = w(i) - l(j) * w(ja(j))
       end do
    end do

    !  Solve U * w = w, where U is upper triangular.
    do i = n, 1, -1
       do j = ua(i) + 1, ia(i+1) - 1
          w(i) = w(i) - l(j) * w(ja(j))
       end do
       w(i) = w(i) / l(ua(i))
    end do

    !  Copy Z out.
    r(1:n) = w(1:n)

    deallocate(w)
  end subroutine lus_cr

  !*****************************************************************************80
  subroutine mult_givens ( c, s, k, g )
!*Description:* This subroutine applies a Givens rotation to two successive entries of a vector.
! In order to make it easier to compare this code with the Original C, the vector indexing is 0-based.
! Input, real(dp) :: C, S, the cosine and sine of a Givens rotation.
!
! Input, integer :: K, indicates the location of the first vector entry.
!
! Input/output, real(dp) :: G(1:K+1), the vector to be modified.
! On output, the Givens rotation has been applied to entries G(K) and G(K+1).

    use arrays, only:dp

    real(dp) :: c
    real(dp) :: s
    integer :: k
    real(dp) :: g(*) !g(1:k+1)

    real(dp) :: g1
    real(dp) :: g2

    g1 = c * g(k) - s * g(k+1)
    g2 = s * g(k) + c * g(k+1)

    g(k)   = g1
    g(k+1) = g2

    return
  end subroutine mult_givens


  !*****************************************************************************80
  subroutine rearrange_cr ( n, nz_num, ia, ja, a )
!*Description:* This subroutine sorts a sparse compressed row matrix.
! It guarantees that the entries in the CR matrix are properly sorted.
!
! After the sorting, the entries of the matrix are rearranged in such
! a way that the entries of each column are listed in ascending order
! of their column values.
! Input, integer :: N, the order of the system.
!
! Input, integer :: NZ_NUM, the number of nonzeros.
!
! Input, integer :: IA(N+1), the compressed row indices.
!
! Input/output, integer :: JA(NZ_NUM), the column indices.
! On output, these may have been rearranged by the sorting.
!
! Input/output, real(dp) :: A(NZ_NUM), the matrix values.  On output,
!    the matrix values may have been moved somewhat because of the sorting.
!
    use arrays, only: dp
    implicit none

    integer :: n
    integer :: nz_num
    integer :: ia(n+1)
    integer :: ja(nz_num)
    real(dp) :: a(nz_num)

    integer :: i
    integer :: itemp
    integer :: k
    integer :: l
    real(dp) :: rtemp

    do i = 1, n

       do k = ia(i), ia(i+1) - 2
          do l = k + 1, ia(i+1) - 1

             if ( ja(l) < ja(k) ) then
                itemp = ja(l)
                ja(l)  = ja(k)
                ja(k)  = itemp

                rtemp = a(l)
                a(l)   = a(k)
                a(k)   = rtemp
             end if

          end do
       end do

    end do

    return
  end subroutine rearrange_cr


end module math_utilities
