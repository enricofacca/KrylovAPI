module IdentityOperators
  use DVectors
  use AbstractOperators
  implicit none
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
    this%nonlinearity=0.d0
    
  end subroutine init_eye

  !! simply copy x into y
  subroutine apply_eye(this,vec_in,vec_out,ierr)
    use DVectors
    implicit none
    class(eye),     intent(inout) :: this
    class(dvector), intent(in   ) :: vec_in
    class(dvector), intent(inout) :: vec_out
    integer,        intent(inout) :: ierr
    
    ! using basic linear algebra operation defined in Dvectors module
    call axpby_dvector(1.0d0,vec_in,0.0d0,vec_out)

    ! or we can simply use a procedure assocaite with standard Fortran vector
    !vec_out%coefficinets=vec_in%coefficients
    
    ierr=0
  end subroutine apply_eye

  
end module IdentityOperators
