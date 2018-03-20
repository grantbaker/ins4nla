module type_defs
  integer, parameter:: sp = kind(1.0),&
    dp = selected_real_kind(2*precision(1.0_sp)),&
    qp = selected_real_kind(2*precision(1.0_dp))
end module type_defs

module problem_setup
  use type_defs
  implicit none
  integer,  parameter :: Nx = 10000
end module problem_setup

module arrs
  use type_defs
  implicit none
  real(dp), allocatable, dimension(:) :: u,b,x
end module arrs

module afuns
  
contains
  
  subroutine apply_1D_laplacian(ax,x,n,hx)
    use type_defs
    implicit none
    integer, intent(in) :: n
    real(dp), intent(out) :: ax(n)
    real(dp), intent(in)  ::  x(n),hx
    real(dp) :: hxi2
    integer :: i
    
    hxi2 = 1.0_dp/hx**2
    
    Ax(1) = hxi2*(x(2) - 2.0_dp*x(1) )
    Ax(n) = hxi2*(     - 2.0_dp*x(n) + x(n-1))
    do i = 2,n-1
     Ax(i)= hxi2*(x(i+1) - 2.0_dp*x(i) + x(i-1))
    end do
    
  end subroutine apply_1D_laplacian
  
end module afuns

module iterative_solvers
  use type_defs
  implicit none
  real(dp), parameter :: TOL = 1.0e-12_dp
  
contains
  
  subroutine jacobi(x,b,nx,hx,l)
    use type_defs
    use afuns
    implicit none
    integer, intent(in) :: nx
    real(dp), intent(inout)  :: hx
    real(dp), intent(inout)  :: x(nx)
    real(dp), intent(inout)  :: b(nx)
    integer, intent(out) :: l
    real(dp) :: ax(nx),xn(nx)
    real(dp) :: delta,rho,s,si
    
    
    delta = TOL*sqrt(sum(b*b))
    rho = 10.0_dp*delta
    si = 0.5_dp*hx**2
    s = 1.0_dp/si
    
    ! loop
    l = 0
    do while (sqrt(rho) .gt. delta)
     call apply_1D_laplacian(ax,x,nx,hx)
     xn = (ax+s*x-b)*si
     rho = sum((x-xn)**2) 
     x = xn
     l = l+1
     ! write(*,*) l, sqrt(rho)
    end do
  end subroutine jacobi
  

subroutine cg(x,b,nx,hx,l)
    use type_defs
    use afuns
    implicit none
    integer, intent(in) :: nx
    real(dp), intent(inout)  :: hx
    real(dp), intent(inout)  :: x(nx)
    real(dp), intent(inout)  :: b(nx)
    integer, intent(out) :: l
    real(dp) :: ax(nx),r(nx),p(nx),q(nx),rtr,alpha,rtrold,beta,A(nx,nx),xn
    real(dp) :: delta,rho,s,si,tol
    integer :: i

    call apply_1D_laplacian(ax,x,nx,hx)
    r=b-ax
    p=r
    rtr=sum(r*r)
!    delta = TOL*sqrt(sum(b*b))
 !   rho = 10.0_dp*delta
 !   si = 0.5_dp*hx**2
 !   s = 1.0_dp/si
    
    ! loop
    l = 0
  !  do while (sqrt(rho) .gt. delta)
    tol=1.0e-5_dp
    do while(sum(r*r) .gt.tol)
       call apply_1D_laplacian(q,p,nx,hx)
       alpha=rtr/sum(p*q)
       x=x+alpha*p
       r=r-alpha*q
       rtrold=rtr
       rtr=sum(r*r)
       beta=rtr/rtrold
       p=r+beta*p
l = l+1
write(*,*) l,rtr
       end do

  end subroutine cg
  
end module iterative_solvers

program ins
  use type_defs
  use problem_setup
  use arrs
  use iterative_solvers
  implicit none
  ! This program solves u_xx = b 
  ! on the domain [x] \in [0,1] with zero boundary conditions 
  ! hx=1/Nx
  real(dp) :: hx
  integer :: i,n_iter  
  ! Set up the grid
  hx = 1.0_dp/real(Nx,dp)
  allocate(x(0:nx))
  do i = 0,nx
   x(i) = real(i,dp)*hx
  end do
  allocate(u(0:nx),b(1:nx-1))
  b = -2.0_dp
  u = 0.0_dp
  call cg(u(1:nx-1),b,nx-1,hx,n_iter)
  write(*,*) "Number of iter ", n_iter, ' Error: '  ,&
    maxval(abs(u - x*(1.0_dp-x)))
end program ins
