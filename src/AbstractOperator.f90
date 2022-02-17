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
  use SimpleVectors
  implicit none
  public :: abstract_operator, application
  public :: eye
  type, abstract :: abstract_operator
     !! Abstract type defining base general operator

     !! Members:     
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
     real(kind=double) :: nonlinearity = 0.0d0
   contains
     !! Procedure for the application of the operator on a vector.
     procedure(application), deferred :: apply
  end type abstract_operator
  
  abstract interface
     !! Abstract procedure defining the interface for a general
     !! application actions
     !!
     !! y = M(x) = ( M*x if M is linear)
     recursive subroutine application(this,vec_in,vec_out,ierr)
       import abstract_operator
       import double
       implicit none
       !! The class with all quantities describing the operator
       class(abstract_operator), intent(inout) :: this
       !! The input" vector 
       real(kind=double),        intent(in   ) :: vec_in(this%ncol)
       !! The output vector y= "operator apllied on" x
       real(kind=double),        intent(inout) :: vec_out(this%nrow)
       !! An integer flag to inform about error
       !! It is zero if no error occured
       integer,                          intent(inout) :: ierr
     end subroutine application
  end interface


  !! identity operator
  type, public, extends(abstract_operator) :: eye
   contains
     !! constructor
     procedure, public, pass :: init => init_eye
     !! procedure to be overrided 
     procedure, public, pass :: apply =>  apply_eye
  end type eye
contains

    !! set identity dimensions
  subroutine init_eye(this,size)
    implicit none
    class(eye),     intent(inout) :: this
    integer,        intent(in   ) :: size
   
    !! the set domain/codomain operator dimension
    this%nrow = size
    this%ncol = size
    this%is_symmetric=.False.
    this%nonlinearity=zero
    
  end subroutine init_eye

  !! simply copy x into y
  subroutine apply_eye(this,vec_in,vec_out,ierr)
    implicit none
    class(eye),        intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: ierr

    ! using basic linear algebra operation defined in SimpleVector module
    call d_axpby(this%nrow,one,vec_in,zero,vec_out)

    ! or we can simply use a procedure assocaite with standard Fortran vector
    vec_out=vec_in
    
    ierr=0
  end subroutine apply_eye

  


  
  
     
end module AbstractOperators
