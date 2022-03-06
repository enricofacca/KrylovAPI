!! Pros:
!!
!!- A user only needs to use the extends(abstract_operator) and define
!!  the apply subroutine to use any Krylov solver.
!!
!!- Defining block matrices, combination of matrices and product
!!  matrices is trivial
!! 
!!- Preconditoners are just some particular abstract operators.
!!
!! Cons.
!! - In same cases, we could save memory and computational power
!! defing other application like y=Ax+y instead of y=Ax. But these
!! savings are small and compromise readibility and flexibility of the
!! approach.
module AbstractOperators
  !! Module containg the abstract class parent of
  !! all operator classes defined on standard Fortran Vector
  use DVectors
  implicit none
  public :: abstract_operator, application
  type, abstract :: abstract_operator
     !! Abstract type defining base general operator

     !! Members shared by any operator:     
     !! Number of columns (the domain of the application)
     integer :: ncol=0
     
     !! Number of rows (the co-domain of the application)
     integer :: nrow=0

     !! Logical flag for symmetry
     logical :: is_symmetric=.false.

     !! Real estimating the non linearity of the application, defined
     !! as: nonlinearity = |OP(x)+OP(y)-OP(x+y)|/(max(|x|,|y|,|x+y|)
     !! It is zero (the defualt value) for linear operators (up to
     !! machine precision). It may be greater that zero for quasi
     !! linear operator like inverse operator.
     double precision :: nonlinearity = 0.0d0
   contains
     !! Procedure for the application of the operator on a vector.
     procedure(application), deferred :: apply
  end type abstract_operator
  
  abstract interface
     !! Abstract procedure defining the interface for a general
     !! application actions
     !!
     !! y = M(x) = ( M*x if M is linear)
     subroutine application(this,vec_in,vec_out,ierr)
       import abstract_operator
       import dvector
       implicit none
       !! The class with all quantities describing the operator
       class(abstract_operator), intent(inout) :: this
       !! The input" vector 
       class(dvector),        intent(in   ) :: vec_in
       !! The output vector y= "operator apllied on" x
       class(dvector),        intent(inout) :: vec_out
       !! An integer flag to inform about error
       !! It is zero if no error occured
       integer,                          intent(inout) :: ierr
     end subroutine application
  end interface    
end module AbstractOperators
