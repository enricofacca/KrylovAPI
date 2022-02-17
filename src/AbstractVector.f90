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

module AbstractVector
  !! Module containg the abstract class parent of
  !! all operator classes.
  !! fixing precision
  implicit none  
  integer, parameter :: double = kind(1.0d0)
  real(kind=double) :: zero=0.0d0
  private
  public :: double, zero
  
  
  type, public, abstract :: abstract_vector
     !! Abstract type defining base general operator

     !! Members:     
     !! Vector size
     integer :: size=0
   contains
     !! Procedure for copy one vector into another
     !procedure(copy), deferred :: copy_procedure
     !! Procedure as daxpy in blas
     !procedure(axpy), deferred :: axpy_procedure
     !! Procedure as dscal, multiplication by a scalar
     !procedure(scal), deferred :: scal_procedure
     !! Procedure scalar product of two vector
     !procedure(scalar_product), deferred :: dot_procedure
     !! Procedure get 
     !procedure(get), deferred :: get_procedure
  end type abstract_vector
  
  ! abstract interface
  !    !! Abstract procedure defining the interface for a general
  !    !! application actions
  !    !!
  !    !! y = M(x) = ( M*x if M is linear)
  !    recursive subroutine application(this,vec_in,vec_out,ierr)
  !      import abstract_operator
  !      import double
  !      implicit none
  !      !! The class with all quantities describing the operator
  !      class(abstract_operator), intent(inout) :: this
  !      !! The input" vector 
  !      real(kind=double),        intent(in   ) :: vec_in(this%ncol)
  !      !! The output vector y= "operator apllied on" x
  !      real(kind=double),        intent(inout) :: vec_out(this%nrow)
  !      !! An integer flag to inform about error
  !      !! It is zero if no error occured
  !      integer,                          intent(inout) :: ierr
  !    end subroutine application
  ! end interface

  
  !! identity operator
  type, public,extends(abstract_vector) :: dvector
     real(kind=double), allocatable :: coefficients(:)
   contains
     !! constructor
     procedure, public, pass :: init => init_dvector
  end type dvector


  !! identity operator
  type, public,extends(abstract_vector) :: dvector_rich
     real(kind=double), allocatable :: coefficients(:)
   contains
     !! constructor
     procedure, public, pass :: init => init_dvector_rich
  end type dvector_rich
contains

  !! set identity dimensions
  subroutine init_dvector(this,size)
    implicit none
    class(dvector),     intent(inout) :: this
    integer,       intent(in   ) :: size
   
    this%size=size
    allocate(this%coefficients(size))
    
  end subroutine init_dvector

  !! set identity dimensions
  subroutine kill_dvector(this,size)
    implicit none
    class(dvector),     intent(inout) :: this
    integer,       intent(in   ) :: size
   
    this%size=size
    deallocate(this%coefficients)
    
  end subroutine kill_dvector

   !! set identity dimensions
  subroutine init_dvector_rich(this,size)
    implicit none
    class(dvector_rich),     intent(inout) :: this
    integer,       intent(in   ) :: size
   
    this%size=size
    allocate(this%coefficients(size))
    
  end subroutine init_dvector_rich

  !! set identity dimensions
  subroutine kill_dvector_rich(this,size)
    implicit none
    class(dvector_rich),     intent(inout) :: this
    integer,       intent(in   ) :: size
   
    this%size=size
    deallocate(this%coefficients)
    
  end subroutine kill_dvector_rich
  
     
end module AbstractVector
