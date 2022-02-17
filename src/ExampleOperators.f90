module ExampleOperators
  use Precision
  use AbstractOperators
  use SimpleVectors

  private

  !! Derived type describing the 1d-Laplacian
  !! with periodic boundary conditions
  type, public, extends(abstract_operator) :: Laplacian
   contains
     procedure, public, pass :: init => init_LP
     !! procedure overrided 
     procedure, public, pass :: apply =>  apply_LP
  end type Laplacian

  !! Derived type describing the 1d-Laplacian
  !! with periodic boundary conditions
  type, public, extends(abstract_operator) :: DiagonalMatrix
     real(kind=double), allocatable :: coefficients(:)
   contains
     procedure, public, pass :: init => init_DM
     procedure, public, pass :: kill => kill_DM
     !! procedure overrided 
     procedure, public, pass :: apply =>  apply_DM
  end type DiagonalMatrix

contains

  


  
  subroutine init_LP(this, nnode)
    class(Laplacian),                intent(inout) :: this
    integer,                          intent(in   ) :: nnode
    
    !! the set domain/codomain operator dimension
    this%nrow = nnode
    this%ncol = nnode
    this%is_symmetric=.True.

  end subroutine init_LP

  subroutine apply_LP(this,vec_in,vec_out,ierr)
    implicit none
    class(Laplacian), intent(inout) :: this
    real(kind=double),  intent(in   ) :: vec_in(this%ncol)
    real(kind=double),  intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: ierr
    !local
    real :: sign=1.0d0
    integer :: i_row


    !call dcopy(this%ncol,vec_in,1,vec_out,1)
    !call dscal(this%ncol,2.0d0,vec_out,1)
    !call daxpy(this%ncol-1,sign,vec_in(2:this%ncol),1,vec_out(1:this%ncol-1),1)
    !call daxpy(this%ncol-1,sign,vec_in(1:this%ncol-1),1,vec_out(2:this%ncol),1)
    
    
    vec_out(1) = 2*vec_in(1) + sign* vec_in(2)
    do i_row=2,this%nrow-1
       vec_out(i_row) = sign* vec_in(i_row-1) + 2*vec_in(i_row) + sign* vec_in(i_row+1)
    end do
    vec_out(this%ncol) =  sign* vec_in(this%ncol-1) + 2*vec_in(this%ncol) 
    ierr=0
  end subroutine apply_LP

  
  !! set dimensions, basic properties and allocate memory
  !! for coefficients on the diagonal
  subroutine init_DM(this, size)
    class(DiagonalMatrix),  intent(inout) :: this
    integer,                intent(in   ) :: size
    
    this%nrow = size
    this%ncol = size
    this%is_symmetric=.True.


    allocate(this%coefficients(size))
    
  end subroutine init_DM

  subroutine kill_DM(this)
    class(DiagonalMatrix),  intent(inout) :: this
    
    this%nrow = 0
    this%ncol = 0
    this%is_symmetric=.False.


    deallocate(this%coefficients)
    
  end subroutine kill_DM

  subroutine apply_DM(this,vec_in,vec_out,ierr)
    implicit none
    class(DiagonalMatrix),  intent(inout) :: this
    real(kind=double), intent(in   ) :: vec_in(this%ncol)
    real(kind=double), intent(inout) :: vec_out(this%nrow)
    integer,           intent(inout) :: ierr
    !local
    integer :: i_row
   
    vec_out=zero
    vec_out=this%coefficients*vec_in
    ierr=0
  end subroutine apply_DM
  
  

end module ExampleOperators
