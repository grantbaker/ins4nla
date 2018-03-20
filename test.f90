module type_defs
  integer, parameter:: sp = kind(1.0),&
    dp = selected_real_kind(2*precision(1.0_sp)),&
    qp = selected_real_kind(2*precision(1.0_dp))
end module type_defs

module new_nla_solvers
contains

  subroutine PCG(A, b)
    use type_defs
    implicit none
    integer :: A, b ! TODO: swap out for real types

    ! TODO: implement preconditioned conjugate gradient method here

  end subroutine PCG

  subroutine precondition_matrix(A, b)
    use type_defs
    implicit none
    integer :: A, b ! TODO: swap out for real types

    ! TODO: precondition A and replace with M^{-1} A

  end subroutine precondition_matrix

  subroutine GS(A, b, n)
    use type_defs
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: A(n,n)
    real(dp), intent(inout) :: b(n)
    real(dp), allocatable, dimension(:) :: x
    integer :: i, j, k
    ! implement Gauss-Seidel and store result in b
    ! TODO: verify that this algorithm is correct

    allocate(x(n))

    do j = 1,n
     x(j) = 0
    end do

    ! TODO: stop iteration at norm, not 10 hardcoded
    do k = 0, 10
     ! set x = b-Ux
     do j = 1,n-1
      x(j) = b(j)
      do i = j+1,n
       x(j) = x(j) - A(i,j)*x(i)
      end do
     x(n) = b(n)
     end do

     ! set x = L^-1 x
     x(1) = x(1)/A(1,1)
     do j = 2, n
      do i = 1, j-1
       x(j) = x(j) - A(i,j)*x(i)
      end do
      x(j) = x(j)/A(j,j)
     end do
    end do

    b = x

  end subroutine GS

end module new_nla_solvers

program test
  use type_defs
  use new_nla_solvers

  real(dp), allocatable, dimension(:) :: x, b
  real(dp), allocatable, dimension(:,:) :: A
  character :: fmt * 8
  integer :: n
  
  n = 3

  allocate(A(n,n))
  allocate(x(n), b(n))

  A = 0
  forall (i = 1:n) A(i,i) = 2
  forall (i = 2:n) A(i,i-1) = 1

  forall (i = 1:n) b(i) = i
  
  write(*, '(A)') 'A'
  write(fmt, '(A,I2,A)') '(', n, 'F6.2)'
  write(*, fmt) A(:,:)
  
  write(*, fmt) b(:)

  call GS(A, b, n)

  write(*, fmt) b(:)

end program test
