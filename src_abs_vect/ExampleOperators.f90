module ExampleOperators
  use Precision
  use AbstractOperators
  use DVectors

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
     double precision, allocatable :: coefficients(:)
   contains
     procedure, public, pass :: init => init_DM
     procedure, public, pass :: kill => kill_DM
     !! procedure overrided 
     procedure, public, pass :: apply =>  apply_DM
  end type DiagonalMatrix
  
contains

  !! constructor
  subroutine init_LP(this, nnode)
    implicit none
    class(Laplacian),                intent(inout) :: this
    integer,                         intent(in   ) :: nnode
    
    !! the set domain/codomain operator dimension
    this%nrow = nnode
    this%ncol = nnode
    this%is_symmetric=.True.

  end subroutine init_LP

  !! apply tridiagonal operator
  subroutine apply_LP(this,vec_in,vec_out,ierr)
    implicit none
    class(Laplacian), intent(inout) :: this
    class(dvector),   intent(in   ) :: vec_in
    class(dvector),   intent(inout) :: vec_out
    integer,          intent(inout) :: ierr
    !local
    double precision :: sign=1.0d0, val
    integer :: i_row

    !! first row
    i_row=1
    val = 2*vec_in%get(i_row) + sign* vec_in%get(i_row+1)
    call vec_out%set(i_row,val)

    !! middle rowsuse DVectors
    do i_row=2,this%nrow-1
       val = sign * vec_in%get(i_row-1) + 2*vec_in%get(i_row) + sign* vec_in%get(i_row+1)
       call vec_out%set(i_row,val)
    end do
    
    !! last row
    i_row=this%nrow
    val = sign* vec_in%get(i_row-1) + 2*vec_in%get(i_row)
    call vec_out%set(i_row,val)

    
    ierr=0
  end subroutine apply_LP
  
  !! constructor
  subroutine init_DM(this, size)
    class(DiagonalMatrix),  intent(inout) :: this
    integer,                intent(in   ) :: size
    
    this%nrow = size
    this%ncol = size
    this%is_symmetric=.True.


    allocate(this%coefficients(size))
    
  end subroutine init_DM

  !! destructor
  subroutine kill_DM(this)
    class(DiagonalMatrix),  intent(inout) :: this
    
    this%nrow = 0
    this%ncol = 0
    this%is_symmetric=.False.

    deallocate(this%coefficients)
    
  end subroutine kill_DM

  !! apply diagonal matrix
  subroutine apply_DM(this,vec_in,vec_out,ierr)
    implicit none
    class(DiagonalMatrix),  intent(inout) :: this
    class(dvector),         intent(in   ) :: vec_in
    class(dvector),         intent(inout) :: vec_out
    integer,                intent(inout) :: ierr
    !local
    integer :: i_row
    double precision :: val

    do i_row=1,this%nrow
       val=vec_in%get(i_row)*this%coefficients(i_row)
       call vec_out%set(i_row,val)
    end do
    ierr=0

  end subroutine apply_DM
  
end module ExampleOperators
