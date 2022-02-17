module SimpleVectors
  !! This module contains the basis linear algebra operation
  !! used in Krylov solver applied to the standard Fortran vectors 
  use Precision
  implicit none
  public

contains
  !! y= alpha*x+beta*y
  !! It should cover daxpy, dcopy, dscal
  subroutine d_axpby(n,alpha,x,beta,y)
    implicit none
    integer,           intent(in   ) :: n
    real(kind=double), intent(in   ) :: alpha
    real(kind=double), intent(in   ) :: x(n)
    real(kind=double), intent(in   ) :: beta
    real(kind=double), intent(inout) :: y(n)

    !! not optimized
    y=alpha*x+beta*y
    
  end subroutine d_axpby

  !! set identity dimensions
  function d_scalar_product(n,x,y) result(xty)
    implicit none
    integer,           intent(in   ) :: n
    real(kind=double), intent(in   ) :: x(n)
    real(kind=double), intent(in   ) :: y(n)
    real(kind=double)                :: xty
    !local
    integer :: i

    !! not optimized
    xty=zero
    do i=1,n
       xty=xty+x(i)*y(i)
    end do
    
  end function d_scalar_product

  !! set identity dimensions
  function d_norm(n,x) result(out)
    implicit none
    integer,           intent(in   ) :: n
    real(kind=double), intent(in   ) :: x(n)
    real(kind=double)                :: out
    out=sqrt(d_scalar_product(n,x,x)) 
  end function d_norm

end module SimpleVectors
