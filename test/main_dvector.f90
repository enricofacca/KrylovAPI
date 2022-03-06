program main
  use AbstractOperators
  !use SimpleVectors
  use DVectors
  use ExampleOperators
 
  use KrylovSolver
  implicit none
  ! problem vars
  type(Laplacian) :: laplacian_1d
  type(DiagonalMatrix) :: inverse_diagonal_laplacian
  type(dvector) :: sol
  type(dvector) :: rhs
  integer :: n
  integer :: max_iterations=400
  double precision :: tolerance=1.0d-5,h,sign

  
  ! local vars
  integer :: ierr,i,m,j
  type(dvector) :: res
  n=100
  h=1.0d0/(n-1)

  !! init matrix
  call laplacian_1d%init(n)

  !! allocate solution and rhs and set the latter
  call sol%init(n)
  call rhs%init(n)
  call res%init(n)

  do i=1,n
     call rhs%set(i,0.d0)
  end do
  call rhs%set(1,1.0d0)
  call rhs%set(n,1.0d0)

  do i=1,n
     call sol%set(i,0.d0)
  end do

  !! init diagonal and set values with inverse of
  !! diagonal of the Laplacian
  call inverse_diagonal_laplacian%init(n)
  inverse_diagonal_laplacian%coefficients=1.0d0/(2.0d0)

  write(*,*) norm_dvector(rhs)
  call pcg(laplacian_1d,rhs,sol,ierr,&
       tolerance=1.0d-3,&
       max_iterations=400)!,&
  !preconditi1.0d0r=inverse_diagonal_laplacian)

  !! check results
 
  call laplacian_1d%apply(sol,res,ierr)
  call axpby_dvector(-1.0d0,rhs,1.0d0,res)
  write(*,*) 'Res=',norm_dvector(res)

end program main
