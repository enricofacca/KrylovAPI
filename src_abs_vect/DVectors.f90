module DVectors
  private
  public :: dvector
  public :: axpby_dvector
  public :: dot_dvector
  public :: norm_dvector

  !! simple real vector 
  type :: dvector
     !! vector length
     integer :: size
     !! vector coefficients
     double precision, allocatable :: coefficients(:)
   contains
     !! constructor
     procedure, public, pass :: init => init_dvector
     !! destructor
     procedure, public, pass :: kill => kill_dvector
     !! axpby procedure 
     procedure, public, nopass :: axpby => axpby_dvector
     !! get procedure 
     procedure, public, pass :: get => get_dvector
     !! set procedure 
     procedure, public, pass :: set => set_dvector
  end type dvector

  !procedure, public :: scalar_product => scalar_product_dvector
  
contains
  !! constructor
  subroutine init_dvector(this,size)
    implicit none
    class(dvector),     intent(inout) :: this
    integer,       intent(in   ) :: size
   
    this%size=size
    allocate(this%coefficients(size))
    
  end subroutine init_dvector

  !! destructor
  subroutine kill_dvector(this)
    implicit none
    class(dvector),     intent(inout) :: this
   
    this%size=0
    deallocate(this%coefficients)
    
  end subroutine kill_dvector

  subroutine axpby_dvector(alpha,x,beta,y)
    use Precision
    implicit none
    double precision,      intent(in   ) :: alpha
    class(dvector),         intent(in   ) :: x
    double precision,      intent(in   ) :: beta
    class(dvector),         intent(inout) :: y

    !! not optimized
    y%coefficients=alpha*x%coefficients+beta*y%coefficients
    
  end subroutine axpby_dvector

  !! get one entry procedure
  function get_dvector(x, index) result(res)
    implicit none
    !! 
    class(dvector),    intent(in   ) :: x
    integer,           intent(in   ) :: index
    double precision  :: res 

    res = x%coefficients(index)
    
  end function get_dvector
  
  !! set one entry procedure
  subroutine set_dvector(x, index, value)
    implicit none
    !! 
    class(dvector),   intent(inout) :: x
    integer,          intent(in   ) :: index
    double precision, intent(in   ) :: value

    x%coefficients(index)=value

  end subroutine set_dvector

  !! scalar product of two vectors
  function dot_dvector(x,y) result(res)
    implicit none
    class(dvector), intent(in   ) :: x
    class(dvector), intent(in   ) :: y
    double precision :: res
    ! local
    integer :: i
    
    res=0.0d0
    do i=1,x%size
       res=res+x%coefficients(i)*y%coefficients(i)
    end do

  end function dot_dvector

  !! vector norm
  function norm_dvector(x) result(res)
    implicit none
    class(dvector), intent(in   ) :: x
    double precision :: res

    res=sqrt(dot_dvector(x,x))
    
  end function norm_dvector

  
end module DVectors
