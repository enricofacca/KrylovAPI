program main
  use AbstractOperator
  use ExampleOperators
  use KrylovSolver
  implicit none
  ! problem vars
  type(Laplacian) :: laplacian_1d
  type(DiagonalMatrix) :: inverse_diagonal_laplacian
  real(kind=double), allocatable :: sol(:)
  real(kind=double), allocatable :: rhs(:)
  integer :: n
  integer :: max_iterations=400
  real(kind=double) :: tolerance=1.0d-5,h,sign
  real(kind=double) :: dnrm2

  
  ! local vars
  integer :: ierr,i,m,j
  real(kind=double), allocatable :: res(:)
  n=100
  h=1.0d0/(n-1)

  !! init matrix
  call laplacian_1d%init(n)

  !! allocate solution and rhs and set the latter
  allocate(sol(n),rhs(n))
  rhs=0
  rhs(1)=1
  rhs(n)=1
  
  sol=0.0d0
  allocate(res(n))
  
  !! init diagonal and set values with inverse of
  !! diagonal of the Laplacian
  call inverse_diagonal_laplacian%init(n)
  inverse_diagonal_laplacian%coefficients=1.0d0/(2.0d0)

  write(*,*) dnrm2(n,rhs,1)
  sol=zero
  call pcg(laplacian_1d,rhs,sol,ierr,&
       tolerance=1.0d-3,&
       max_iterations=400)!,&
  !preconditioner=inverse_diagonal_laplacian)

  !! check results
 
  call laplacian_1d%apply(sol,res,ierr)
  res=res-rhs;
  write(*,*) 'Res=',dnrm2(n,res,1)

end program main
