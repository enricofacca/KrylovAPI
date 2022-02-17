module TransposableOperators
  use AbstractOperators
  !! Abstract type defining base for linear operator for which
  !! the transpose operator che be easily defined.

  !! make calss and method public
  public :: transposable_operator, application_transpose
  
  type, abstract, extends(abstract_operator) :: transposable_operator
     !! No further elements are require.
     !! Non lienarity is zero (by default)
   contains
     !! Procedure for computation of the operator applied on a vector
     !! which the a matrix-vector multiplication for linear operators.
     procedure(application_transpose), deferred :: apply_transpose
  end type transposable_operator

  abstract interface
     !! Abstract procedure defining the interface for a general
     !! application actions
     !!
     !! y = M^T x
     subroutine application_transpose(this,x,y,ierr)
       import transposable_operator
       import double
       implicit none
       class(transposable_operator), target, intent(inout) :: this
       real(kind=double), target,            intent(in   ) :: x(this%nrow)
       real(kind=double), target,            intent(inout) :: y(this%ncol)
       integer,                              intent(inout) :: ierr
     end subroutine application_transpose
  end interface
end module TransposableOperators
